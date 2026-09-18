"""Driver: replaces the Fortran ``regionalised_daily`` entry point.

Parameters arrive as arguments rather than through ``data_r.dat``, and the
generated series is returned as an array rather than written to
``mmm_<idx>.out`` and read back.  ``nearby_station_details.out`` is still read
from disk, because :func:`pyraingen.getnearbystations.station` writes it and
users edit it by hand.
"""

from datetime import date

import numpy as np

from . import preprocess, simulate
from .calendarutil import NDAY, NOUT, daycount
from .io import read_nearby_stations, read_station_data
from .poolstats import build_pools
from .rng import new_state, warm_up

WARMUP_CALLS = 500
DEFAULT_SEED = 131


def _cumulative_weights(weight, nstn):
    """Cumulative donor-station weights, with the last forced to exactly 1."""
    ws = np.zeros(nstn + 1, dtype=np.float32)
    for i in range(1, nstn + 1):
        ws[i] = weight[i]
        if i > 1:
            ws[i] = ws[i] + ws[i - 1]
    ws[nstn] = np.float32(1.0)
    return ws


def _simulated_days(startyear, nyears):
    """Number of days in the simulated period, honouring leap years."""
    total = 0
    nsg = startyear - 1
    for i in range(1, nyears + 1):
        for j in range(1, NOUT + 1):
            nd = int(NDAY[j])
            if j == 2:
                nd = daycount(nsg, i)
            total += nd
    return total


def run(data_path, nearby_path, nyears, startyear, nsim, nstation=5,
        cutoff=0.30, wind=15, nlon=3, lag=1, iamt=1, z=None, nout=NOUT,
        seed=DEFAULT_SEED, verbose=True):
    """Generate daily rainfall for the target site.

    Parameters
    ----------
    data_path : str
        Directory holding the ``rev_drNNNNN.txt`` daily files.
    nearby_path : str
        Path to ``nearby_station_details.out``.
    nyears, startyear, nsim : int
        Length, start year and number of realisations.
    nstation : int
        Number of nearby stations to borrow statistics from.
    cutoff : float
        Wet-day threshold in mm.
    wind : int
        Half-width of the day window, in days.
    nlon : int
        Number of long-memory windows.
    lag : int
        Order of the occurrence Markov chain.  Only ``1`` is supported, which
        is the only value the Fortran was ever run with.
    iamt : int
        1 to generate amounts, 0 for wet/dry flags only.
    z : sequence
        ``[[starts], [ends]]`` for the long-memory windows, in days.
    seed : int
        Seed for ``ran1``.  Note that the generator reseeds with
        ``max(-seed, 1)``, so every value ``>= 0`` gives the same stream as the
        historical default of 131; pass a negative value for a different one.

    Returns
    -------
    tuple
        ``(rain, day_vector)`` with ``rain`` shaped ``(ndays, nsim)`` in mm and
        ``day_vector`` the matching Julian day numbers.
    """
    if z is None:
        z = [[2, 2, 2], [90, 180, 345]]
    if lag != 1:
        raise NotImplementedError(
            "only lag=1 is supported; the Fortran was never run with lag=0"
        )

    rain = np.float32(cutoff)
    iband = int(wind)
    nstn = int(nstation)
    nlon = int(nlon)

    lon = np.zeros((3, nlon + 1), dtype=np.int64)
    for k in range(nlon):
        lon[1, k + 1] = int(z[0][k])
        lon[2, k + 1] = int(z[1][k])

    if nlon == 0:
        iyrst = 1
    else:
        longest = int(lon[2, 1:].max())
        iyrst = int(longest / 365.0 + 2.001)

    state, gset = new_state(seed)
    warm_up(state, gset, WARMUP_CALLS)

    if verbose:
        print("Reading data")
    station_ids, weight = read_nearby_stations(nearby_path, nstn)
    rf, nyrs, nstrt, _avrf = read_station_data(data_path, station_ids, nstn, nout)
    ws = _cumulative_weights(weight, nstn)

    nyear_max = int(nyrs[1:].max())
    nbuf = (2 * iband + 1) * nyear_max + 2
    nday = NDAY

    if verbose:
        print("Calculating probabilities .....")
    pro1, pro2, pro = preprocess.smoothprob(rf, nday, nstn, nyrs, nstrt, iband,
                                            rain, nout, nbuf)

    if nlon > 0:
        high = preprocess.long_store(rf, nday, nstn, nyrs, nstrt, nlon, lon,
                                     iyrst, rain, nout, nyear_max)
        avl, sdl = preprocess.av_sd_lon(rf, high, nday, nstn, nyrs, nstrt,
                                        nlon, iyrst, iband, nout, nbuf)
        if verbose:
            print("Calculating means and variances of long term variables")
        lav, lcov, plmin, plmax = preprocess.smoothavcov(
            rf, high, sdl, pro1, pro2, nday, nstn, nyrs, nstrt, nlon, iyrst,
            iband, rain, nout, nbuf)
        del high
    else:
        sdl = np.ones((nlon + 1, 13, 32, nstn + 1), dtype=np.float32)
        lav = np.zeros((2, 2, 13, 32, nlon + 2, nstn + 1), dtype=np.float32)
        lcov = np.zeros((2, 2, 13, 32, nlon * nlon + 1, nstn + 1),
                        dtype=np.float32)
        plmin = np.zeros((13, 32, nlon + 1, nstn + 1), dtype=np.float32)
        plmax = np.zeros((13, 32, nlon + 1, nstn + 1), dtype=np.float32)

    if verbose:
        print("Storing variables for rf_amount")
    counts = preprocess.amt_counts(rf, nday, nstn, nyrs, nstrt, iband, rain,
                                   nout)
    offsets = (np.cumsum(counts.ravel()) - counts.ravel()).reshape(
        counts.shape).astype(np.int64)
    total = int(counts.sum())
    pool_y = np.zeros(max(total, 1), dtype=np.float32)
    pool_x = np.zeros(max(total, 1), dtype=np.float32)
    pool_next = np.zeros(max(total, 1), dtype=np.int64)
    preprocess.amt_fill(rf, nday, nstn, nyrs, nstrt, iband, rain, nout,
                        offsets, counts, pool_y, pool_x, pool_next, nbuf)

    (cls_y, cls_x, coffsets, ccounts, hgy, hgx, sxxi_t, sxxixy_t, scond_t,
     sd_t, ymax_t) = build_pools(pool_y, pool_x, pool_next, offsets, counts,
                                 rain, nday, nstn, nout)
    del pool_y, pool_x, pool_next

    ndays = _simulated_days(startyear, nyears)
    maxpool = int(ccounts.max())
    nsg = startyear - 1

    if verbose:
        print("Simulating rainfall.......")
    out = simulate.run(state, gset, rf, nday, nstn, nyrs, nstrt, ws, nlon, lon,
                       iyrst, rain, nout, int(nsim), int(nyears), nsg, pro1,
                       pro2, lav, lcov, plmin, plmax, sdl, cls_y, cls_x,
                       coffsets, ccounts, hgy, hgx, sxxi_t, sxxixy_t, scond_t,
                       sd_t, ymax_t, ndays, maxpool)

    if iamt == 0:
        out = (out >= rain).astype(np.float32)

    day_vector = _day_vector(startyear, ndays)
    return out.T.copy(), day_vector


def _day_vector(startyear, ndays):
    """Julian day numbers for the simulated period."""
    from ..datevectojd import datevecToJD

    return np.arange(0, ndays, 1) + datevecToJD(date(int(startyear), 1, 1))
