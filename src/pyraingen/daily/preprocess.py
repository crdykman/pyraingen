"""Preprocessing: the once-per-run tables the day loop reads from.

Ports ``smoothprob``, ``long_store``/``prevdays``, ``av_sd_lon``,
``smoothavcov``/``prevavsd1``/``solve`` and ``rf_amt_store``/``rank1``.

Every kernel here is a deliberate transliteration rather than a tidy-up.  The
Fortran routines disagree with one another in ways that change results: the
leap-year epoch ``ns`` is ``nstrt-1`` in ``smoothprob``, ``smoothavcov`` and
``long_store`` but ``nstrt`` in ``av_sd_lon`` and ``rf_amt_store``, and the
year loop starts at 1 in some and at ``iyrst`` in others.  Sums are
accumulated sequentially in float32 for the same reason: pairwise summation
would round differently and the difference propagates through the chain.
"""

import math

import numpy as np

from ._jit import njit
from .calendarutil import day_neg, day_pos, daycount, prev_day

F0 = np.float32(0.0)
F1 = np.float32(1.0)


@njit(cache=True)
def basic(data, n):
    """Mean and standard deviation, matching the Fortran ``BASIC``."""
    ave = F0
    sd = F1
    if n <= 1:
        return ave, sd
    for j in range(n):
        ave += data[j]
    ave = ave / np.float32(n)
    var = F0
    ep = F0
    for j in range(n):
        s = data[j] - ave
        ep += s
        var += s * s
    sd = np.float32(math.sqrt((var - ep * ep / np.float32(n)) / np.float32(n - 1)))
    return ave, sd


# ---------------------------------------------------------------------------
# smoothprob
# ---------------------------------------------------------------------------


@njit(cache=True)
def smoothprob(rf, nday, nstn, nyrs, nstrt, iband, rain, nout, nbuf):
    """Window-pooled wet-day probabilities for every (month, day, station).

    Returns ``(pro1, pro2, pro)`` -- P(wet | previous dry), P(wet | previous
    wet) and P(wet) -- each shaped ``(13, 32, nstn+1)``.
    """
    pro1 = np.zeros((13, 32, nstn + 1), dtype=np.float32)
    pro2 = np.zeros((13, 32, nstn + 1), dtype=np.float32)
    pro = np.zeros((13, 32, nstn + 1), dtype=np.float32)

    xc = np.zeros(nbuf, dtype=np.float32)
    xp1 = np.zeros(nbuf, dtype=np.float32)

    for j in range(1, nout + 1):
        for l in range(1, nday[j] + 1):
            for nk in range(1, nstn + 1):
                kk = 0
                nyear = nyrs[nk]
                ns = nstrt[nk] - 1

                for i in range(1, nyear + 1):
                    if j == 2 and l > daycount(ns, i):
                        continue

                    li = l - iband - 1
                    for _jh in range(2 * iband + 1):
                        li += 1
                        ic = i
                        jc = j
                        lc = li
                        if li <= 2:
                            ic, jc, lc, indx = day_neg(ic, jc, lc, 1, ns, 2, nday)
                        else:
                            ic, jc, lc, indx = day_pos(ic, jc, lc, ns, nyear, nday)
                        if indx == 1:
                            continue

                        ip, jp, lp, bad = prev_day(ic, jc, lc, ns, nday)
                        if bad == 1:
                            continue

                        # Drop the entry when tomorrow falls past end of record.
                        i_n = ic
                        jn = jc
                        ln = lc + 1
                        nd = nday[jn]
                        if jn == 2:
                            nd = daycount(ns, i_n)
                        if i_n == nyear and jn == nout and ln > nd:
                            continue
                        if ln > nd:
                            jn += 1
                            ln = 1
                            if jn > nout:
                                i_n += 1
                                jn = 1
                            if i_n > nyear:
                                continue

                        xc[kk] = rf[nk, ic, jc, lc]
                        xp1[kk] = rf[nk, ip, jp, lp]
                        kk += 1

                # prob1
                sum1 = F0
                sum2 = F0
                sum5 = F0
                k1 = 0
                k2 = 0
                for i in range(kk):
                    if xc[i] >= rain and xp1[i] < rain:
                        sum1 += F1
                    if xc[i] >= rain and xp1[i] >= rain:
                        sum2 += F1
                    if xc[i] >= rain:
                        sum5 += F1
                    if xp1[i] < rain:
                        k1 += 1
                    else:
                        k2 += 1
                if k1 > 0:
                    pro1[j, l, nk] = sum1 / np.float32(k1)
                if k2 > 0:
                    pro2[j, l, nk] = sum2 / np.float32(k2)
                if kk > 0:
                    pro[j, l, nk] = sum5 / np.float32(kk)

    return pro1, pro2, pro


# ---------------------------------------------------------------------------
# long_store / prevdays
# ---------------------------------------------------------------------------


@njit(cache=True)
def prevdays(rf, nday, i1, j1, l1, long1, long2, nk, ns, rain, nout):
    """Count wet days between ``long1`` and ``long2`` days before (i1, j1, l1).

    Returns ``(total, ind)``; ``ind == 1`` means the walk ran off the start of
    the record.
    """
    total = F0
    lp = l1
    jp = j1
    ip = i1
    for l in range(1, long2 + 1):
        lp -= 1
        if lp < 1:
            if ip == 1 and jp == 1:
                return total, 1
            jp -= 1
            if jp < 1:
                ip -= 1
                if ip < 1:
                    return total, 1
                jp = nout
            lp = nday[jp]
            if jp == 2:
                lp = daycount(ns, ip)
        if l >= long1 and rf[nk, ip, jp, lp] >= rain:
            total += F1
    return total, 0


@njit(cache=True)
def long_store(rf, nday, nstn, nyrs, nstrt, nlon, lon, iyrst, rain, nout,
               nyear_max):
    """Wet-day counts in each long-memory window, for every day of record.

    Returns ``high`` shaped ``(nstn*nlon+1, nyear_max+2, 13, 32)``.

    The two-stage initialisation matters and is not tidy-up-able.  In the
    Fortran ``high`` lives in a static common block, so it starts as *zero*,
    and the routine's own init loop writes the -9999 "no value" marker only
    over days that exist in the calendar.  29 February of a non-leap year is
    therefore left at 0.0, and the ``high < 0`` filter in ``av_sd_lon`` and
    ``smoothavcov`` lets it straight through: those days enter the pools as a
    spurious zero long-memory count.  Filling the whole array with -9999
    instead would silently drop them and shift every February statistic.
    """
    high = np.zeros((nstn * nlon + 1, nyear_max + 2, 13, 32), dtype=np.float32)

    for nk in range(1, nstn + 1):
        nyear = nyrs[nk]
        ns = nstrt[nk] - 1
        for i in range(1, nyear + 1):
            for j in range(1, nout + 1):
                lc = nday[j]
                if j == 2:
                    lc = daycount(ns, i)
                for l in range(1, lc + 1):
                    for jj in range(1, nlon + 1):
                        jk = (nk - 1) * nlon + jj
                        high[jk, i, j, l] = np.float32(-9999.0)

    for nk in range(1, nstn + 1):
        nyear = nyrs[nk]
        ns = nstrt[nk] - 1
        for i in range(iyrst, nyear + 1):
            for j in range(1, nout + 1):
                lc = nday[j]
                if j == 2:
                    lc = daycount(ns, i)
                for l in range(1, lc + 1):
                    for jj in range(1, nlon + 1):
                        dsum, ind = prevdays(rf, nday, i, j, l, lon[1, jj],
                                             lon[2, jj], nk, ns, rain, nout)
                        if ind == 1:
                            break
                        jk = (nk - 1) * nlon + jj
                        high[jk, i, j, l] = dsum

    return high


# ---------------------------------------------------------------------------
# av_sd_lon
# ---------------------------------------------------------------------------


@njit(cache=True)
def av_sd_lon(rf, high, nday, nstn, nyrs, nstrt, nlon, iyrst, iband, nout,
              nbuf):
    """Window-pooled mean and SD of each long-memory count.

    Returns ``(avl, sdl)`` shaped ``(nlon+1, 13, 32, nstn+1)``.
    """
    avl = np.zeros((nlon + 1, 13, 32, nstn + 1), dtype=np.float32)
    sdl = np.zeros((nlon + 1, 13, 32, nstn + 1), dtype=np.float32)

    x1 = np.zeros((nbuf, nlon + 1), dtype=np.float32)
    col = np.zeros(nbuf, dtype=np.float32)

    for j in range(1, nout + 1):
        for l in range(1, nday[j] + 1):
            for nk in range(1, nstn + 1):
                nn = 0
                nyear = nyrs[nk]
                ns = nstrt[nk]

                for i in range(iyrst, nyear + 1):
                    if j == 2 and l > daycount(ns, i):
                        continue

                    li = l - iband - 1
                    for _jh in range(2 * iband + 1):
                        li += 1
                        ic = i
                        jc = j
                        lc = li
                        if li <= 2:
                            ic, jc, lc, indx = day_neg(ic, jc, lc, 1, ns, 2, nday)
                        else:
                            ic, jc, lc, indx = day_pos(ic, jc, lc, ns, nyear, nday)
                        if indx == 1:
                            continue

                        _ip, _jp, _lp, bad = prev_day(ic, jc, lc, ns, nday)
                        if bad == 1:
                            continue

                        skip = False
                        for jj in range(1, nlon + 1):
                            jk = (nk - 1) * nlon + jj
                            if high[jk, ic, jc, lc] < F0:
                                skip = True
                                break
                        if skip:
                            continue

                        for jj in range(1, nlon + 1):
                            jk = (nk - 1) * nlon + jj
                            x1[nn, jj] = high[jk, ic, jc, lc]
                        nn += 1

                for jj in range(1, nlon + 1):
                    for i in range(nn):
                        col[i] = x1[i, jj]
                    ave, sd = basic(col, nn)
                    avl[jj, j, l, nk] = ave
                    sdl[jj, j, l, nk] = sd
                    if j == 2 and l == 29:
                        avl[jj, j, l, nk] = avl[jj, j, l - 1, nk]
                        sdl[jj, j, l, nk] = sdl[jj, j, l - 1, nk]

    return avl, sdl


# ---------------------------------------------------------------------------
# solve (Numerical Recipes SVD pseudo-inverse)
# ---------------------------------------------------------------------------


@njit(cache=True)
def _jacobi_eigh(a, nv):
    """Cyclic Jacobi eigendecomposition of a symmetric matrix.

    Returns ``(eigenvalues, eigenvectors)`` with eigenvectors in columns.
    Rotations are done in float64 for stability; the caller casts back.
    """
    v = np.zeros((nv, nv), dtype=np.float64)
    for i in range(nv):
        v[i, i] = 1.0

    for _sweep in range(60):
        off = 0.0
        for p in range(nv - 1):
            for q in range(p + 1, nv):
                off += a[p, q] * a[p, q]
        if off <= 1.0e-30:
            break

        for p in range(nv - 1):
            for q in range(p + 1, nv):
                if abs(a[p, q]) <= 1.0e-300:
                    continue
                theta = (a[q, q] - a[p, p]) / (2.0 * a[p, q])
                if theta >= 0.0:
                    t = 1.0 / (theta + math.sqrt(theta * theta + 1.0))
                else:
                    t = -1.0 / (-theta + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c

                for k in range(nv):
                    akp = a[k, p]
                    akq = a[k, q]
                    a[k, p] = c * akp - s * akq
                    a[k, q] = s * akp + c * akq
                for k in range(nv):
                    apk = a[p, k]
                    aqk = a[q, k]
                    a[p, k] = c * apk - s * aqk
                    a[q, k] = s * apk + c * aqk
                for k in range(nv):
                    vkp = v[k, p]
                    vkq = v[k, q]
                    v[k, p] = c * vkp - s * vkq
                    v[k, q] = s * vkp + c * vkq

    w = np.zeros(nv, dtype=np.float64)
    for i in range(nv):
        w[i] = a[i, i]
    return w, v


@njit(cache=True)
def solve(ss, nv):
    """Overwrite ``ss`` with its pseudo-inverse, returning the determinant.

    Mirrors the Fortran ``solve``: an *absolute* singular-value cutoff of
    1e-7, and a determinant taken as the product of the singular values, so
    always non-negative.  ``nv == 1`` takes the reciprocal directly and leaves
    a zero variance as zero rather than inverting it -- both behaviours matter
    downstream.

    The original went through Numerical Recipes' ``svdcmp``.  Every matrix
    reaching here is a covariance matrix and therefore symmetric, so a Jacobi
    eigendecomposition gives the same decomposition: for a symmetric ``A`` with
    eigenpairs ``(lambda_i, q_i)`` the SVD is ``w_i = |lambda_i|``, ``V = Q``
    and ``U = Q diag(sign lambda)``, which makes ``V diag(1/w) U^T`` exactly
    ``Q diag(1/lambda) Q^T``.  Doing it this way keeps the kernel free of any
    LAPACK dependency, which numba's ``np.linalg`` would otherwise drag in.
    """
    tol = np.float32(1.0e-7)

    if nv == 1:
        det = ss[0, 0]
        if ss[0, 0] != F0:
            ss[0, 0] = F1 / ss[0, 0]
        return det

    a = np.zeros((nv, nv), dtype=np.float64)
    for i in range(nv):
        for j in range(nv):
            a[i, j] = np.float64(ss[i, j])

    lam, q = _jacobi_eigh(a, nv)

    for i in range(nv):
        for j in range(nv):
            acc = 0.0
            for k in range(nv):
                if abs(lam[k]) > tol:
                    acc += q[i, k] * q[j, k] / lam[k]
            ss[i, j] = np.float32(acc)

    det = F1
    for i in range(nv):
        det = det * np.float32(abs(lam[i]))
    return det


# ---------------------------------------------------------------------------
# smoothavcov / prevavsd1
# ---------------------------------------------------------------------------


@njit(cache=True)
def prevavsd1(xc, xp1, xpl, kp, nv, rain):
    """Class-conditional means and pooled inverse covariance of the
    standardised long-memory sums.

    Splits the window pool by (today wet?, yesterday wet?) and returns
    ``(xbar, sighat_inv, det)`` where ``xbar`` is indexed
    ``[today, yesterday, variable]`` with 0-based 0=dry, 1=wet.

    Called with ``iflag=1`` by :func:`smoothavcov`, so all four classes share
    the pooled covariance ``sighat``; only the means differ.
    """
    xbar = np.zeros((2, 2, nv), dtype=np.float32)
    xss = np.zeros((2, 2, nv, nv), dtype=np.float32)
    nnn = np.zeros((2, 2), dtype=np.int64)

    for i in range(kp):
        ixc = 1 if xc[i] >= rain else 0
        ixp = 1 if xp1[i] >= rain else 0
        nnn[ixc, ixp] += 1
        for k in range(nv):
            xbar[ixc, ixp, k] += xpl[i, k + 1]
            for ll in range(nv):
                xss[ixc, ixp, k, ll] += xpl[i, k + 1] * xpl[i, ll + 1]

    for i in range(2):
        for jj in range(2):
            if nnn[i, jj] > 0:
                for k in range(nv):
                    xbar[i, jj, k] = xbar[i, jj, k] / np.float32(nnn[i, jj])

    df = np.float32(kp - 4)

    sighat = np.zeros((nv, nv), dtype=np.float32)
    for k in range(nv):
        for ll in range(nv):
            acc = F0
            for i in range(2):
                for jj in range(2):
                    acc += (xss[i, jj, k, ll]
                            - np.float32(nnn[i, jj]) * xbar[i, jj, k]
                            * xbar[i, jj, ll])
            if df > F0:
                acc = acc / df
            sighat[k, ll] = acc

    det = solve(sighat, nv)
    if kp < 10:
        det = F0

    return xbar, sighat, det


@njit(cache=True)
def smoothavcov(rf, high, sdl, pro1, pro2, nday, nstn, nyrs, nstrt, nlon,
                iyrst, iband, rain, nout, nbuf):
    """Means, inverse covariances and pool ranges of the long-memory sums.

    Returns ``(lav, lcov, plmin, plmax)``:

    * ``lav[today, yesterday, month, day, 1..nlon, station]`` -- class means,
      with the determinant parked in slot ``nlon+1`` exactly as the Fortran
      does.
    * ``lcov[today, yesterday, month, day, flat 3x3, station]``.
    * ``plmin``/``plmax`` -- the pooled range of each window count.

    ``today``/``yesterday`` are 0 for dry and 1 for wet.  February 29 copies
    February 28, and ``pro1``/``pro2`` are patched in place for that day as
    the Fortran does at the end of ``smoothavcov``.
    """
    lav = np.zeros((2, 2, 13, 32, nlon + 2, nstn + 1), dtype=np.float32)
    lcov = np.zeros((2, 2, 13, 32, nlon * nlon + 1, nstn + 1), dtype=np.float32)
    plmin = np.zeros((13, 32, nlon + 1, nstn + 1), dtype=np.float32)
    plmax = np.zeros((13, 32, nlon + 1, nstn + 1), dtype=np.float32)

    xc = np.zeros(nbuf, dtype=np.float32)
    xp1 = np.zeros(nbuf, dtype=np.float32)
    xpl = np.zeros((nbuf, nlon + 1), dtype=np.float32)

    for j in range(1, nout + 1):
        for l in range(1, nday[j] + 1):
            for nk in range(1, nstn + 1):
                pmin = np.full(nlon + 1, np.float32(10000.0), dtype=np.float32)
                pmax = np.full(nlon + 1, np.float32(-10000.0), dtype=np.float32)

                kk = 0
                nyear = nyrs[nk]
                ns = nstrt[nk] - 1

                for i in range(iyrst, nyear + 1):
                    if j == 2 and l > daycount(ns, i):
                        continue

                    li = l - iband - 1
                    for _jh in range(2 * iband + 1):
                        li += 1
                        ic = i
                        jc = j
                        lc = li
                        if li <= 2:
                            ic, jc, lc, indx = day_neg(ic, jc, lc, 1, ns, 2, nday)
                        else:
                            ic, jc, lc, indx = day_pos(ic, jc, lc, ns, nyear, nday)
                        if indx == 1:
                            continue

                        ip, jp, lp, bad = prev_day(ic, jc, lc, ns, nday)
                        if bad == 1:
                            continue

                        skip = False
                        for jj in range(1, nlon + 1):
                            jk = (nk - 1) * nlon + jj
                            if high[jk, ic, jc, lc] < F0:
                                skip = True
                                break
                        if skip:
                            continue

                        for jj in range(1, nlon + 1):
                            jk = (nk - 1) * nlon + jj
                            aa = high[jk, ic, jc, lc]
                            if aa < pmin[jj]:
                                pmin[jj] = aa
                            if aa > pmax[jj]:
                                pmax[jj] = aa
                            xpl[kk, jj] = aa / sdl[jj, j, l, nk]

                        xc[kk] = rf[nk, ic, jc, lc]
                        xp1[kk] = rf[nk, ip, jp, lp]
                        kk += 1

                xbar, sighat, det = prevavsd1(xc, xp1, xpl, kk, nlon, rain)

                for icur in range(2):
                    for iprev in range(2):
                        j3 = 0
                        for j1 in range(nlon):
                            for j2 in range(nlon):
                                j3 += 1
                                lcov[icur, iprev, j, l, j3, nk] = sighat[j1, j2]
                        for j1 in range(nlon):
                            lav[icur, iprev, j, l, j1 + 1, nk] = xbar[icur, iprev, j1]
                        lav[icur, iprev, j, l, nlon + 1, nk] = det

                for j1 in range(1, nlon + 1):
                    plmin[j, l, j1, nk] = pmin[j1]
                    plmax[j, l, j1, nk] = pmax[j1]

    # February 29 copies February 28 throughout, and the occurrence
    # probabilities are averaged with 1 March.
    for nk in range(1, nstn + 1):
        for icur in range(2):
            for iprev in range(2):
                for j3 in range(1, nlon * nlon + 1):
                    lcov[icur, iprev, 2, 29, j3, nk] = lcov[icur, iprev, 2, 28, j3, nk]
                for j3 in range(1, nlon + 2):
                    lav[icur, iprev, 2, 29, j3, nk] = lav[icur, iprev, 2, 28, j3, nk]
        for j1 in range(1, nlon + 1):
            plmin[2, 29, j1, nk] = plmin[2, 28, j1, nk]
            plmax[2, 29, j1, nk] = plmax[2, 28, j1, nk]
        pro1[2, 29, nk] = (pro1[2, 28, nk] + pro1[3, 1, nk]) / np.float32(2.0)
        pro2[2, 29, nk] = (pro2[2, 28, nk] + pro2[3, 1, nk]) / np.float32(2.0)

    return lav, lcov, plmin, plmax


# ---------------------------------------------------------------------------
# rf_amt_store / rank1
# ---------------------------------------------------------------------------


@njit(cache=True)
def rank1(xc, xp, ixn, n):
    """Sort the pool into descending order of ``xc``, carrying ``xp``/``ixn``.

    This is the Fortran ``rank1`` bubble/selection sort transliterated rather
    than replaced by ``argsort``.  It is *not* stable -- with a pool of
    0.1 mm-resolution rainfall depths there are many ties, and which tied
    analogue ends up where decides both the paired previous-day amount and the
    inverse-CDF pick.  Using a stable sort here changes the generated series.
    """
    for i in range(n - 1):
        for j in range(i + 1, n):
            if xc[j] > xc[i]:
                a = xc[i]
                xc[i] = xc[j]
                xc[j] = a
                a = xp[i]
                xp[i] = xp[j]
                xp[j] = a
                ia = ixn[i]
                ixn[i] = ixn[j]
                ixn[j] = ia


@njit(cache=True)
def amt_counts(rf, nday, nstn, nyrs, nstrt, iband, rain, nout):
    """Number of wet analogues in each (month, day, station) window pool."""
    counts = np.zeros((13, 32, nstn + 1), dtype=np.int64)

    for j in range(1, nout + 1):
        for l in range(1, nday[j] + 1):
            for nk in range(1, nstn + 1):
                nn = 0
                nyear = nyrs[nk]
                ns = nstrt[nk]

                for i in range(1, nyear + 1):
                    if j == 2 and l > daycount(ns, i):
                        continue

                    li = l - iband - 1
                    for _jh in range(2 * iband + 1):
                        li += 1
                        ic = i
                        jc = j
                        lc = li
                        if li <= 2:
                            ic, jc, lc, indx = day_neg(ic, jc, lc, 1, ns, 2, nday)
                        else:
                            ic, jc, lc, indx = day_pos(ic, jc, lc, ns, nyear, nday)
                        if indx == 1:
                            continue

                        _ip, _jp, _lp, bad = prev_day(ic, jc, lc, ns, nday)
                        if bad == 1:
                            continue

                        i_n = ic
                        jn = jc
                        ln = lc + 1
                        nd = nday[jn]
                        if jn == 2:
                            nd = daycount(ns, i_n)
                        if i_n == nyear and jn == nout and ln > nd:
                            continue
                        if ln > nd:
                            jn += 1
                            ln = 1
                            if jn > nout:
                                i_n += 1
                                jn = 1
                            if i_n > nyear:
                                continue

                        if rf[nk, ic, jc, lc] >= rain:
                            nn += 1

                counts[j, l, nk] = nn

    return counts


@njit(cache=True)
def amt_fill(rf, nday, nstn, nyrs, nstrt, iband, rain, nout, offsets, counts,
             pool_y, pool_x, pool_next, nbuf):
    """Fill the flat, descending-sorted analogue pools.

    ``pool_y`` holds today's wet-day amount, ``pool_x`` yesterday's amount and
    ``pool_next`` a flag for whether tomorrow was wet.
    """
    xc = np.zeros(nbuf, dtype=np.float32)
    xp = np.zeros(nbuf, dtype=np.float32)
    ixn = np.zeros(nbuf, dtype=np.int64)

    for j in range(1, nout + 1):
        for l in range(1, nday[j] + 1):
            for nk in range(1, nstn + 1):
                nn = 0
                nyear = nyrs[nk]
                ns = nstrt[nk]

                for i in range(1, nyear + 1):
                    if j == 2 and l > daycount(ns, i):
                        continue

                    li = l - iband - 1
                    for _jh in range(2 * iband + 1):
                        li += 1
                        ic = i
                        jc = j
                        lc = li
                        if li <= 2:
                            ic, jc, lc, indx = day_neg(ic, jc, lc, 1, ns, 2, nday)
                        else:
                            ic, jc, lc, indx = day_pos(ic, jc, lc, ns, nyear, nday)
                        if indx == 1:
                            continue

                        ip, jp, lp, bad = prev_day(ic, jc, lc, ns, nday)
                        if bad == 1:
                            continue

                        i_n = ic
                        jn = jc
                        ln = lc + 1
                        nd = nday[jn]
                        if jn == 2:
                            nd = daycount(ns, i_n)
                        if i_n == nyear and jn == nout and ln > nd:
                            continue
                        if ln > nd:
                            jn += 1
                            ln = 1
                            if jn > nout:
                                i_n += 1
                                jn = 1
                            if i_n > nyear:
                                continue

                        irhn = 0
                        if rf[nk, i_n, jn, ln] >= rain:
                            irhn = 1

                        if rf[nk, ic, jc, lc] >= rain:
                            xc[nn] = rf[nk, ic, jc, lc]
                            xp[nn] = rf[nk, ip, jp, lp]
                            ixn[nn] = irhn
                            nn += 1

                rank1(xc, xp, ixn, nn)

                base = offsets[j, l, nk]
                for ii in range(nn):
                    pool_y[base + ii] = xc[ii]
                    pool_x[base + ii] = xp[ii]
                    pool_next[base + ii] = ixn[ii]
