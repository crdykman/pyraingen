"""Tests for the IFD conditioning path.

These cover the defects fixed in 2.1.0. The whole module was untestable before
that release: ``import pyraingen.ifdcond`` raised ``AttributeError`` on
numpy >= 2 because ``np.RankWarning`` had been removed, so nothing downstream
of the import had ever been exercised on a modern numpy.

The end-to-end tests build their own small sub-daily netCDF rather than using
the 32 MB bundled example. That keeps them quick, and because the synthetic
input is sampled hourly (24 records per day) rather than every six minutes, it
also exercises the fix that made ``computeIFD`` and ``aggregateRainfall``
honour ``minsPerSample`` instead of hardcoding 240 records and 6 minutes.
"""

from datetime import date

import numpy as np
import pytest

from pyraingen.aggregaterainfall import aggregateRainfall
from pyraingen.computeifd import computeIFD
from pyraingen.datevectojd import datevecToJD
from pyraingen.ifdcond import ifdcond
from pyraingen.ifdcondobjfun import ifdcondobjfun
from pyraingen.producesubdailynetcdf import produceSubDailyNetCDF

# Hourly sampling keeps the synthetic files small and exercises the
# non-default sampling path.
RECORDS_PER_DAY = 24
MINS_PER_SAMPLE = 60
YEAR_START, YEAR_END = 2001, 2012
N_SIMS_IN_FILE = 5
AEP = [63.20, 50, 20, 10]        # four points for the cubic polyfit
DURATIONS = [60, 120]            # 1 and 2 records at hourly sampling


def _n_days(year_start, year_end):
    return (date(year_end, 12, 31) - date(year_start, 1, 1)).days + 1


@pytest.fixture(scope="module")
def synthetic_input(tmp_path_factory):
    """A small sub-daily netCDF in the schema readSynthRainNetCDF expects."""
    out = tmp_path_factory.mktemp("ifd") / "subdaily.nc"
    n_days = _n_days(YEAR_START, YEAR_END)
    rng = np.random.default_rng(11)

    # Intermittent rainfall: mostly dry, with occasional wet hours.
    data = rng.random((RECORDS_PER_DAY, n_days, N_SIMS_IN_FILE))
    data[data < 0.85] = 0.0
    data *= 12.0

    day_vector = np.arange(0, n_days) + datevecToJD(date(YEAR_START, 1, 1))
    produceSubDailyNetCDF(str(out), data.astype(np.float32), day_vector,
                          title="synthetic", institution="test")
    return str(out)


@pytest.fixture(scope="module")
def target_ifd(tmp_path_factory):
    """Target IFD: rows are AEP, columns are durations."""
    path = tmp_path_factory.mktemp("ifd_target") / "target.csv"
    rows = [[18.0, 26.0], [15.0, 22.0], [11.0, 16.0], [8.0, 12.0]]
    path.write_text("\n".join(",".join(f"{v}" for v in r) for r in rows) + "\n")
    return str(path)


@pytest.fixture(scope="module")
def single_duration_target(tmp_path_factory):
    """A one-column target IFD, which numpy loads as 1-D without ndmin=2."""
    path = tmp_path_factory.mktemp("ifd_one") / "target_one.csv"
    path.write_text("18.0\n15.0\n11.0\n8.0\n")
    return str(path)


# ---------------------------------------------------------------------------
# aggregateRainfall
# ---------------------------------------------------------------------------


def test_aggregation_sums_consecutive_records():
    rng = np.random.default_rng(1)
    series = rng.random((3, 20, 240))
    got = aggregateRainfall(series, 10)

    assert got.shape == (3, 20, 24)
    assert np.allclose(got, series.reshape(3, 20, 24, 10).sum(axis=3))


def test_aggregation_derives_records_per_day_from_the_data():
    """It used to hardcode 240, silently mis-aggregating other sampling rates."""
    series = np.random.default_rng(2).random((2, 5, RECORDS_PER_DAY))
    got = aggregateRainfall(series, 2)

    assert got.shape == (2, 5, RECORDS_PER_DAY // 2)
    assert np.allclose(got, series.reshape(2, 5, 12, 2).sum(axis=3))


def test_aggregation_refuses_an_uneven_division():
    """It used to truncate, silently discarding the remainder records."""
    series = np.random.default_rng(3).random((2, 4, 240))
    with pytest.raises(ValueError, match="does not divide evenly"):
        aggregateRainfall(series, 7)


def test_aggregation_conserves_total_depth():
    series = np.random.default_rng(4).random((2, 6, 240))
    assert aggregateRainfall(series, 10).sum() == pytest.approx(series.sum())


# ---------------------------------------------------------------------------
# computeIFD
# ---------------------------------------------------------------------------


def test_computeifd_extracts_annual_maxima():
    n_years, days_per_year = 3, 10
    years = np.repeat([2001, 2002, 2003], days_per_year)
    series = np.zeros((2, n_years * days_per_year, RECORDS_PER_DAY))
    # Plant a known maximum in each year of each simulation.
    planted = {}
    for sim in range(2):
        for i, year in enumerate([2001, 2002, 2003]):
            value = 10.0 + sim + i
            series[sim, i * days_per_year + 2, 5] = value
            planted[(sim, year)] = value

    ifd = computeIFD(series, years, np.array([MINS_PER_SAMPLE]),
                     minsPerSample=MINS_PER_SAMPLE)

    assert ifd.shape == (n_years, 2, 1)
    for sim in range(2):
        expected = sorted(planted[(sim, y)] for y in [2001, 2002, 2003])
        assert np.allclose(ifd[:, sim, 0], expected)


def test_computeifd_sorts_ascending_over_years():
    rng = np.random.default_rng(6)
    years = np.repeat([2001, 2002, 2003, 2004], 12)
    series = rng.random((3, 48, RECORDS_PER_DAY))
    ifd = computeIFD(series, years, np.array([MINS_PER_SAMPLE, 2 * MINS_PER_SAMPLE]),
                     minsPerSample=MINS_PER_SAMPLE)

    assert ifd.shape == (4, 3, 2)
    assert np.all(np.diff(ifd, axis=0) >= 0)


def test_computeifd_honours_minspersample():
    """With the wrong sampling assumption the aggregation window is wrong."""
    rng = np.random.default_rng(7)
    years = np.repeat([2001, 2002], 15)
    series = rng.random((2, 30, RECORDS_PER_DAY))

    hourly = computeIFD(series, years, np.array([120]),
                        minsPerSample=MINS_PER_SAMPLE)          # 2 records
    manual = aggregateRainfall(series, 2)
    expected = np.array([
        [manual[s, years == y, :].max() for s in range(2)]
        for y in [2001, 2002]
    ])
    assert np.allclose(np.sort(expected, axis=0), hourly[:, :, 0])


# ---------------------------------------------------------------------------
# ifdcondobjfun
# ---------------------------------------------------------------------------


def test_objective_is_zero_for_a_perfect_match():
    n_years, n_sims, n_dur = 20, 3, 2
    rng = np.random.default_rng(8)
    simulated = np.sort(rng.random((n_years, n_sims, n_dur)) * 20, axis=0)
    freq = np.array([63, 50, 20, 10])
    # Build a target equal to the mean across simulations at the matched ranks.
    idx = np.array([np.argmin(abs(f - 100 * np.arange(1, n_years + 1) / (n_years + 1)))
                    for f in freq])
    target = np.stack([simulated[idx, :, d].mean(axis=1) for d in range(n_dur)],
                      axis=1)

    rmae = ifdcondobjfun(target, simulated, freq)
    assert rmae.shape == (n_dur,)
    assert np.allclose(rmae, 0.0, atol=1e-12)


def test_objective_handles_a_single_simulation():
    """np.squeeze used to collapse the simulation axis, raising an AxisError."""
    n_years, n_dur = 20, 2
    rng = np.random.default_rng(9)
    simulated = np.sort(rng.random((n_years, 1, n_dur)) * 20, axis=0)
    target = np.full((4, n_dur), 10.0)

    rmae = ifdcondobjfun(target, simulated, np.array([63, 50, 20, 10]))
    assert rmae.shape == (n_dur,)
    assert np.isfinite(rmae).all()


# ---------------------------------------------------------------------------
# ifdcond argument validation (fails before any heavy computation)
# ---------------------------------------------------------------------------


def test_duration_absent_from_the_estimated_set_is_rejected(synthetic_input,
                                                            target_ifd,
                                                            tmp_path):
    with pytest.raises(ValueError, match="are not present in"):
        ifdcond(synthetic_input, str(tmp_path / "out.nc"), target_ifd,
                nSims=2, nRecursions=1,
                TargetIFDdurationsEst=DURATIONS,
                TargetIFDdurations=[60, 999],
                AEP=AEP, nRecordsPerDay=RECORDS_PER_DAY,
                minsPerSample=MINS_PER_SAMPLE, plot=False, seed=1)


def test_a_proper_subset_of_durations_is_rejected(synthetic_input, target_ifd,
                                                  tmp_path):
    """The docstring used to promise subsets worked; they raised IndexError.

    The limitation is now stated and enforced up front instead.
    """
    with pytest.raises(NotImplementedError, match="must currently equal"):
        ifdcond(synthetic_input, str(tmp_path / "out.nc"), target_ifd,
                nSims=2, nRecursions=1,
                TargetIFDdurationsEst=DURATIONS,
                TargetIFDdurations=[60],
                AEP=AEP, nRecordsPerDay=RECORDS_PER_DAY,
                minsPerSample=MINS_PER_SAMPLE, plot=False, seed=1)


# ---------------------------------------------------------------------------
# ifdcond end to end
# ---------------------------------------------------------------------------


def _read_rainfall(path):
    import netCDF4 as nc

    ds = nc.Dataset(path)
    data = ds["rainfall"][:].data
    ds.close()
    return data


def test_ifdcond_runs_and_writes_valid_output(synthetic_input, target_ifd,
                                              tmp_path):
    out = tmp_path / "conditioned.nc"
    ifdcond(synthetic_input, str(out), target_ifd,
            nSims=2, nRecursions=1,
            TargetIFDdurationsEst=DURATIONS, TargetIFDdurations=DURATIONS,
            AEP=AEP, nRecordsPerDay=RECORDS_PER_DAY,
            minsPerSample=MINS_PER_SAMPLE, plot=False, seed=1)

    assert out.exists()
    data = _read_rainfall(str(out))
    assert data.shape[0] == RECORDS_PER_DAY
    assert data.shape[2] == 2
    assert np.isfinite(data).all()
    assert (data >= 0).all()
    assert data.sum() > 0


def test_ifdcond_is_reproducible_for_a_given_seed(synthetic_input, target_ifd,
                                                  tmp_path):
    """The simulation subsample used the global RNG, so runs never repeated."""
    outputs = []
    for run in range(2):
        out = tmp_path / f"seeded_{run}.nc"
        ifdcond(synthetic_input, str(out), target_ifd,
                nSims=2, nRecursions=1,
                TargetIFDdurationsEst=DURATIONS, TargetIFDdurations=DURATIONS,
                AEP=AEP, nRecordsPerDay=RECORDS_PER_DAY,
                minsPerSample=MINS_PER_SAMPLE, plot=False, seed=1234)
        outputs.append(_read_rainfall(str(out)))

    assert np.array_equal(outputs[0], outputs[1])


def test_a_single_duration_run_completes(synthetic_input,
                                         single_duration_target, tmp_path):
    """Used to crash twice over.

    ``plt.subplots(nrows=1)`` returns a bare Axes so ``ax[i]`` raised, and
    ``np.loadtxt`` returns a 1-D array for a one-column target IFD so every
    ``targetIFD[:, column]`` access failed. Plotting is left on here because
    the subplots bug is only reachable with plot=True.
    """
    out = tmp_path / "one_duration.nc"
    cwd = tmp_path / "plots"
    cwd.mkdir()

    import os
    previous = os.getcwd()
    os.chdir(cwd)                      # the figure is written to the cwd
    try:
        ifdcond(synthetic_input, str(out), single_duration_target,
                nSims=2, nRecursions=1,
                TargetIFDdurationsEst=[60], TargetIFDdurations=[60],
                AEP=AEP, nRecordsPerDay=RECORDS_PER_DAY,
                minsPerSample=MINS_PER_SAMPLE, plot=True, seed=2)
    finally:
        os.chdir(previous)

    assert out.exists()
    assert (cwd / "ifdcond_plot.png").exists()
    data = _read_rainfall(str(out))
    assert np.isfinite(data).all()
    assert (data >= 0).all()
