"""Analogue pools and the statistics that depend only on them.

The Fortran recomputed ``rank_h``, ``hfracx``, ``hfracy`` and the whole
covariance block inside ``psimmain`` on *every* wet day -- together roughly
70% of its runtime -- even though all of it depends only on the pool for
(day, station, class) and not on the conditioning value or the random draws.

Here those quantities are computed once per pool.  The per-wet-day work then
drops to O(n): form the weights, cumulate them, pick one analogue.

Pools are stored flat with an offset/length table so they stay numba-friendly.
The class index follows the Fortran ``icg``/``ih`` convention:

====  ==============  =============
icg   yesterday       tomorrow
====  ==============  =============
1     wet             wet
2     wet             dry
3     dry             wet
4     dry             dry
====  ==============  =============

Classes 3 and 4 always take the unconditional branch: with a dry previous day
there is no conditioning variable, so ``nvv`` is 0 and ``ind`` is forced to 1.
"""

import math

import numpy as np

from ._jit import njit

F0 = np.float32(0.0)
F1 = np.float32(1.0)
HUPPER = np.float32(3.5)
HLOWER = np.float32(1.0) / np.float32(3.5)

# Minimum pool size for the conditional KDE; below this ``rf_amt_gen`` sets
# ind=1 and falls back to unconditional resampling.
NMIN_COND = 15


@njit(cache=True)
def rank_h(v, n):
    """Interquartile-ish range used as the kernel window width.

    Ports the Fortran ``rank_h``: sort descending, then take the gap between
    the 25th and 75th order statistics at the original's 1-based indices
    ``int(n*0.25+0.05)`` and ``int(n*0.75+0.05)``.  Those can be 0 for tiny
    pools, which would index out of bounds in Fortran, so they are clamped.
    """
    x = np.sort(v[:n])
    i1 = int(np.float32(n) * np.float32(0.25) + np.float32(0.05))
    i2 = int(np.float32(n) * np.float32(0.75) + np.float32(0.05))
    if i1 < 1:
        i1 = 1
    if i2 < 1:
        i2 = 1
    if i1 > n:
        i1 = n
    if i2 > n:
        i2 = n
    # x is ascending; the Fortran indexed a descending sort.
    return abs(x[n - i1] - x[n - i2])


@njit(cache=True)
def hfrac(v, n, out):
    """Per-point bandwidth factors, ports of ``hfracx``/``hfracy``.

    For each point a window of +/- 0.5 * ``rank_h`` is placed around it, slid
    back inside the pool range if it overhangs, and the number of pool members
    inside decides the bandwidth.

    The window is clipped against ``amx``/``amn``, which the Fortran never
    assigned -- they were stack garbage, so bandwidths near the pool extremes
    were compiler- and call-history-dependent.  Here they are the pool maximum
    and minimum, the evident intent.  This changes results relative to
    pyraingen 1.0.2 regardless of the port.

    Note that ``href`` here omits the 0.8 factor that ``psimmain`` applies to
    its own ``h``; that asymmetry is in the original and is kept.
    """
    nv = 1
    fac = np.float32((4.0 / (nv + 2.0)) ** (1.0 / (nv + 4.0)))
    href = fac * np.float32(np.float32(n) ** np.float32(-1.0 / (nv + 4.0)))

    sd = rank_h(v, n)

    amx = v[0]
    amn = v[0]
    for i in range(n):
        if v[i] > amx:
            amx = v[i]
        if v[i] < amn:
            amn = v[i]

    for i in range(n):
        aup = v[i] + np.float32(0.5) * sd
        alr = v[i] - np.float32(0.5) * sd
        if aup > amx:
            alr = alr - (aup - amx)
            aup = amx
        if alr < amn:
            aup = aup + (amn - alr)
            alr = amn

        ii = 0
        for jj in range(n):
            if v[jj] > alr and v[jj] < aup:
                ii += 1

        g = (np.float32(1.2) - np.float32(ii) / np.float32(n)) * href
        if g > HUPPER * href:
            g = HUPPER * href
        if g < HLOWER * href:
            g = HLOWER * href
        out[i] = g


@njit(cache=True)
def split_counts(pool_x, pool_next, offsets, counts, rain, nday, nstn, nout):
    """Size of each of the four analogue classes per (month, day, station)."""
    ccounts = np.zeros((13, 32, nstn + 1, 5), dtype=np.int64)
    for j in range(1, nout + 1):
        for l in range(1, nday[j] + 1):
            for nk in range(1, nstn + 1):
                base = offsets[j, l, nk]
                for i in range(counts[j, l, nk]):
                    if pool_x[base + i] >= rain:
                        ih = 1 if pool_next[base + i] > 0 else 2
                    else:
                        ih = 3 if pool_next[base + i] > 0 else 4
                    ccounts[j, l, nk, ih] += 1
    return ccounts


@njit(cache=True)
def split_fill(pool_y, pool_x, pool_next, offsets, counts, rain, nday, nstn,
               nout, coffsets, cls_y, cls_x):
    """Split each pool into its four classes, preserving descending order."""
    cursor = np.zeros((13, 32, nstn + 1, 5), dtype=np.int64)
    for j in range(1, nout + 1):
        for l in range(1, nday[j] + 1):
            for nk in range(1, nstn + 1):
                base = offsets[j, l, nk]
                for i in range(counts[j, l, nk]):
                    if pool_x[base + i] >= rain:
                        ih = 1 if pool_next[base + i] > 0 else 2
                    else:
                        ih = 3 if pool_next[base + i] > 0 else 4
                    at = coffsets[j, l, nk, ih] + cursor[j, l, nk, ih]
                    cls_y[at] = pool_y[base + i]
                    cls_x[at] = pool_x[base + i]
                    cursor[j, l, nk, ih] += 1


@njit(cache=True)
def pool_statistics(cls_y, cls_x, coffsets, ccounts, nday, nstn, nout,
                    hgy, hgx, sxxi_t, sxxixy_t, scond_t, sd_t, ymax_t):
    """Precompute everything in ``psimmain`` that does not depend on the day.

    Fills, per (month, day, station, class):

    * ``hgy`` -- ``hfracy`` bandwidth factors for every candidate analogue.
    * ``hgx`` -- ``hfracx`` bandwidth factors (conditional classes only).
    * ``sxxi_t``/``sxxixy_t``/``scond_t`` -- the regression of today's amount
      on yesterday's, from ``estxcpdfstats`` and ``estscond``.
    * ``sd_t`` -- pool SD, used by the unconditional branch.
    * ``ymax_t`` -- the largest analogue, used by the 3x cap.
    """
    for j in range(1, nout + 1):
        for l in range(1, nday[j] + 1):
            for nk in range(1, nstn + 1):
                for ih in range(1, 5):
                    n = ccounts[j, l, nk, ih]
                    if n <= 0:
                        continue
                    base = coffsets[j, l, nk, ih]
                    y = cls_y[base:base + n]
                    x = cls_x[base:base + n]

                    ymax_t[j, l, nk, ih] = y[0]

                    # basic(y, av, sd, n) for the unconditional branch
                    ave = F0
                    if n > 1:
                        for i in range(n):
                            ave += y[i]
                        ave = ave / np.float32(n)
                        var = F0
                        ep = F0
                        for i in range(n):
                            s = y[i] - ave
                            ep += s
                            var += s * s
                        sd_t[j, l, nk, ih] = np.float32(
                            math.sqrt((var - ep * ep / np.float32(n))
                                      / np.float32(n - 1)))
                    else:
                        sd_t[j, l, nk, ih] = F1

                    hfrac(y, n, hgy[base:base + n])

                    # Classes 3 and 4 have no conditioning variable, and short
                    # pools fall back to unconditional resampling, so the
                    # regression block is never reached for them.
                    if ih > 2 or n < NMIN_COND:
                        continue

                    hfrac(x, n, hgx[base:base + n])

                    xmn = F0
                    for i in range(n):
                        xmn += x[i]
                    xmn = xmn / np.float32(n)

                    cov = F0
                    for i in range(n):
                        d = x[i] - xmn
                        cov += d * d
                    cov = cov / np.float32(n)

                    # solve() with nv == 1: a zero variance is left as zero.
                    sxxi = cov
                    if cov != F0:
                        sxxi = F1 / cov
                    sxxi_t[j, l, nk, ih] = sxxi

                    ymn = F0
                    for i in range(n):
                        ymn += y[i]
                    ymn = ymn / np.float32(n)

                    syy = F0
                    sxy = F0
                    for i in range(n):
                        dy = y[i] - ymn
                        syy += dy * dy
                        sxy += (x[i] - xmn) * dy
                    syy = syy / np.float32(n)
                    sxy = sxy / np.float32(n)

                    sxxixy = sxxi * sxy
                    scond = syy - sxy * sxxixy
                    if scond < F0:
                        scond = F0
                    sxxixy_t[j, l, nk, ih] = sxxixy
                    scond_t[j, l, nk, ih] = scond


def build_pools(pool_y, pool_x, pool_next, offsets, counts, rain, nday, nstn,
                nout):
    """Split the sorted pools by class and precompute their statistics.

    Returns a tuple of arrays consumed by :func:`pyraingen.daily.simulate.run`.
    """
    ccounts = split_counts(pool_x, pool_next, offsets, counts, rain, nday,
                           nstn, nout)

    flat = np.cumsum(ccounts.ravel()) - ccounts.ravel()
    coffsets = flat.reshape(ccounts.shape).astype(np.int64)
    total = int(ccounts.sum())

    cls_y = np.zeros(max(total, 1), dtype=np.float32)
    cls_x = np.zeros(max(total, 1), dtype=np.float32)
    split_fill(pool_y, pool_x, pool_next, offsets, counts, rain, nday, nstn,
               nout, coffsets, cls_y, cls_x)

    hgy = np.zeros(max(total, 1), dtype=np.float32)
    hgx = np.zeros(max(total, 1), dtype=np.float32)
    shape = (13, 32, nstn + 1, 5)
    sxxi_t = np.zeros(shape, dtype=np.float32)
    sxxixy_t = np.zeros(shape, dtype=np.float32)
    scond_t = np.zeros(shape, dtype=np.float32)
    sd_t = np.zeros(shape, dtype=np.float32)
    ymax_t = np.zeros(shape, dtype=np.float32)

    pool_statistics(cls_y, cls_x, coffsets, ccounts, nday, nstn, nout, hgy,
                    hgx, sxxi_t, sxxixy_t, scond_t, sd_t, ymax_t)

    return (cls_y, cls_x, coffsets, ccounts, hgy, hgx, sxxi_t, sxxixy_t,
            scond_t, sd_t, ymax_t)
