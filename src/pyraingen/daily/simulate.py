"""The serial day loop.

Nothing here can be vectorised.  There is a single ``ran1`` stream driving
station choice, occurrence and amounts; realisations do not restart but
continue the Markov state of the one before; and each day's amount conditions
on the previous day's generated depth.  So the loop stays a loop, and the
ordering of random draws is part of the specification.

Two quirks of the original are reproduced deliberately and are called out
where they occur: the one-day lag in what gets written, and the fact that the
multivariate-normal densities are never normalised.
"""

import math

import numpy as np

from ._jit import njit
from .calendarutil import daycount
from .poolstats import NMIN_COND
from .rng import gasdev, ran1

F0 = np.float32(0.0)
F1 = np.float32(1.0)
RAIN_CAP = np.float32(600.0)


@njit(cache=True)
def firsttimeh(rf, nday, ist, nk, nstrt, nlon, lon, rain, nout, sumlon):
    """Seed the long-memory ring buffers from a historical year.

    Fills ``sumlon[k, 1..lon[2,k]]`` with the wet flags of the ``lon[2,k]``
    days ending on 31 December of year slot ``ist``.  Returns 1 if the record
    does not reach back far enough.
    """
    ns = nstrt[nk] - 1
    for kk in range(1, nlon + 1):
        ln = lon[2, kk]
        im = ist
        jm = nout
        nd = nday[jm]
        if jm == 2:
            nd = daycount(ns, im)
        ls = nd

        while True:
            js = ls - ln
            if js > 0:
                break
            jm -= 1
            if jm < 1:
                jm = nout
                im -= 1
            if im < 1:
                return 1
            nd = nday[jm]
            if jm == 2:
                nd = daycount(ns, im)
            ls = ls + nd

        js -= 1
        nd = nday[jm]
        if jm == 2:
            nd = daycount(ns, im)

        for l in range(1, ln + 1):
            js += 1
            if js > nd:
                js -= nd
                jm += 1
                if jm > nout:
                    jm = 1
                    im += 1
                nd = nday[jm]
                if jm == 2:
                    nd = daycount(ns, im)
            sumlon[kk, l] = F0
            if rf[nk, im, jm, js] >= rain:
                sumlon[kk, l] = F1
    return 0


@njit(cache=True)
def check_lon(igen, sumlon, lon, nlon):
    """Shift the ring buffers on by a day and refresh the window sums.

    Yesterday's occurrence flag is pushed in at the newest slot, and the sum
    over each window is parked in slot ``lon[2,k]+1``.
    """
    for kk in range(1, nlon + 1):
        nn = lon[2, kk]
        if nn > 1:
            for l in range(1, nn):
                sumlon[kk, l] = sumlon[kk, l + 1]
        sumlon[kk, nn] = np.float32(igen)

    for kk in range(1, nlon + 1):
        nn = lon[2, kk] - lon[1, kk] + 1
        total = F0
        for l in range(1, nn + 1):
            total += sumlon[kk, l]
        sumlon[kk, lon[2, kk] + 1] = total


@njit(cache=True)
def multiply(x, av, cov, nv):
    """Unnormalised multivariate normal kernel, accumulated in float64."""
    qf = 0.0
    for k in range(nv):
        for l in range(nv):
            qf += (float(cov[k, l]) * (float(x[k]) - float(av[k]))
                   * (float(x[l]) - float(av[l])))
    part = -0.5 * qf
    if part < -15.0:
        part = -15.0
    return np.float32(math.exp(part))


@njit(cache=True)
def rf_amt_gen(state, gset, jc, lc, nk, p1, inext, rain, cls_y, cls_x,
               coffsets, ccounts, hgy, hgx, sxxi_t, sxxixy_t, scond_t, sd_t,
               ymax_t, wts):
    """Draw one wet-day amount.

    ``vv`` and ``xr`` are consumed *before* the pool is formed, and exactly one
    of each is consumed on every wet day whichever branch is taken, so the
    stream stays in lockstep with the Fortran.
    """
    vv = ran1(state)
    xr = gasdev(state, gset)

    if p1 >= rain:
        icg = 1 if inext > 0 else 2
    else:
        icg = 3 if inext > 0 else 4

    n = ccounts[jc, lc, nk, icg]
    if n <= 0:
        return rain

    base = coffsets[jc, lc, nk, icg]
    ymax = ymax_t[jc, lc, nk, icg]

    # Classes 3 and 4 carry no conditioning variable, so nvv is 0 there.
    nvv = 0 if icg > 2 else 1
    unconditional = nvv < 1 or n < NMIN_COND

    if unconditional:
        iwt = int(vv * np.float32(n) + F1)
        if iwt < 1:
            iwt = 1
        if iwt > n:
            iwt = n
        if n < 5:
            ys = cls_y[base + iwt - 1]
        else:
            h = hgy[base + iwt - 1]
            ys = cls_y[base + iwt - 1] + h * sd_t[jc, lc, nk, icg] * xr
            if ys > ymax * np.float32(3.0):
                ys = cls_y[base + iwt - 1]
    else:
        sxxi = sxxi_t[jc, lc, nk, icg]
        sumwt = F0
        for i in range(n):
            xx = p1 - cls_x[base + i]
            hg = hgx[base + i]
            uu = xx * sxxi / (hg * hg) * xx
            uu1 = uu / np.float32(2.0)
            if uu1 > np.float32(20.0):
                uu1 = F0
            else:
                uu1 = np.float32(math.exp(-uu1))
            wts[i] = uu1
            sumwt += uu1

        if sumwt > F0:
            for i in range(n):
                wts[i] = wts[i] / sumwt
        else:
            for i in range(n):
                wts[i] = F1 / np.float32(n)

        for i in range(1, n):
            wts[i] = wts[i - 1] + wts[i]

        rnum = vv * wts[n - 1]
        iwt = n
        for i in range(n):
            if wts[i] >= rnum:
                iwt = i + 1
                break

        scond = scond_t[jc, lc, nk, icg]
        hused = hgy[base + iwt - 1]
        b = ((p1 - cls_x[base + iwt - 1]) * sxxixy_t[jc, lc, nk, icg]
             * hused / hgx[base + iwt - 1])
        bran = hused * np.float32(math.sqrt(scond)) * xr
        ys = cls_y[base + iwt - 1] + b + bran

        if ys > ymax * np.float32(3.0):
            if b > np.float32(1.5) * ymax:
                b = np.float32(1.5) * ymax
            ys = cls_y[base + iwt - 1] + b + bran

        ii = 0
        while True:
            ii += 1
            if ii > 5:
                break
            if ys < F0:
                hused = hused / np.float32(2.0)
                bran = hused * np.float32(math.sqrt(scond)) * xr
                if b < F0:
                    b = b / np.float32(2.0)
                ys = cls_y[base + iwt - 1] + b + bran
            else:
                break

    if ys > RAIN_CAP:
        ys = RAIN_CAP
    if ys < rain:
        ys = rain
    return ys


@njit(cache=True)
def run(state, gset, rf, nday, nstn, nyrs, nstrt, ws, nlon, lon, iyrst, rain,
        nout, nsim, ng, nsg, pro1, pro2, lav, lcov, plmin, plmax, sdl,
        cls_y, cls_x, coffsets, ccounts, hgy, hgx, sxxi_t, sxxixy_t, scond_t,
        sd_t, ymax_t, ndays, maxpool):
    """Generate ``nsim`` realisations of ``ng`` years.

    Returns ``out`` shaped ``(nsim, ndays)`` in chronological order.
    """
    out = np.zeros((nsim, ndays), dtype=np.float32)

    maxlon = 0
    for kk in range(1, nlon + 1):
        if lon[2, kk] > maxlon:
            maxlon = lon[2, kk]
    sumlon = np.zeros((nlon + 1, maxlon + 2), dtype=np.float32)

    wts = np.zeros(max(maxpool, 1), dtype=np.float32)
    prevlong = np.zeros(nlon, dtype=np.float32)
    pav = np.zeros(nlon, dtype=np.float32)
    ppav = np.zeros(nlon, dtype=np.float32)
    pcov = np.zeros((nlon, nlon), dtype=np.float32)
    ppcov = np.zeros((nlon, nlon), dtype=np.float32)

    # --- initial state: a donor station and a historical year to seed from ---
    vv = ran1(state)
    nk = nstn
    for kk in range(1, nstn + 1):
        if vv <= ws[kk]:
            nk = kk
            break
    nyear = nyrs[nk]

    while True:
        ist = int(ran1(state) * np.float32(nyear) + F1)
        if ist < iyrst or ist > nyear:
            continue
        ind = 0
        if nlon > 0:
            ind = firsttimeh(rf, nday, ist, nk, nstrt, nlon, lon, rain, nout,
                             sumlon)
        if ind == 0:
            break

    nd = nday[nout]
    ip1 = F0
    if rf[nk, ist, nout, nd] >= rain:
        ip1 = F1

    igen = 0
    if rf[nk, ist, nout, 31] >= rain:
        igen = 1
    p1 = rf[nk, ist, nout, 31]
    gen = rf[nk, ist, nout, 31]
    icur = igen

    jp = 1
    lp = 1

    for it in range(nsim):
        w = 0
        for i in range(1, ng + 1):
            for j in range(1, nout + 1):
                nd = nday[j]
                if j == 2:
                    nd = daycount(nsg, i)

                for l1 in range(1, nd + 1):
                    # donor station for today's statistics
                    vv = ran1(state)
                    nk = nstn
                    for kk in range(1, nstn + 1):
                        if vv <= ws[kk]:
                            nk = kk
                            break

                    if nlon > 0:
                        check_lon(igen, sumlon, lon, nlon)

                    if ip1 < rain:
                        pw = pro1[j, l1, nk]
                        iprev = 0
                    else:
                        pw = pro2[j, l1, nk]
                        iprev = 1

                    pp = pw
                    if nlon > 0:
                        for j1 in range(nlon):
                            aa = sumlon[j1 + 1, lon[2, j1 + 1] + 1]
                            prevlong[j1] = aa / sdl[j1 + 1, j, l1, nk]

                        detl = lav[0, iprev, j, l1, nlon + 1, nk]
                        ddetl = lav[1, iprev, j, l1, nlon + 1, nk]

                        if detl != F0 and ddetl != F0:
                            j3 = 0
                            for j1 in range(nlon):
                                for j2 in range(nlon):
                                    j3 += 1
                                    pcov[j1, j2] = lcov[0, iprev, j, l1, j3, nk]
                                    ppcov[j1, j2] = lcov[1, iprev, j, l1, j3, nk]
                                pav[j1] = lav[0, iprev, j, l1, j1 + 1, nk]
                                ppav[j1] = lav[1, iprev, j, l1, j1 + 1, nk]

                            pd = multiply(prevlong, pav, pcov, nlon)
                            pn = multiply(prevlong, ppav, ppcov, nlon)

                            # The Fortran guards the (2*pi)^(-nlon/2)*det^(-1/2)
                            # normalisation with `det`/`ddet`, which are never
                            # assigned -- so it never fires.  Reproduced: the
                            # factors do not cancel unless det == ddet.

                            for jj in range(nlon):
                                aa = sumlon[jj + 1, lon[2, jj + 1] + 1]
                                if (aa < (F1 / np.float32(1.5))
                                        * plmin[j, l1, jj + 1, nk]
                                        or aa > np.float32(1.5)
                                        * plmax[j, l1, jj + 1, nk]):
                                    pd = F1
                                    pn = F1

                            pp = pn * pp / (pp * pn + (F1 - pp) * pd)

                    pw = pp

                    # occurrence
                    vv = ran1(state)
                    igen = 0
                    if vv <= pw:
                        igen = 1

                    # The amount written now is for the *previous* day: the
                    # kernel conditions on whether the next day is wet, which
                    # is only known once today's occurrence has been drawn.
                    # Only the month and day are needed, to look up that day's
                    # analogue pool; the year is implicit in the write order.
                    jc = jp
                    lc = lp
                    inext = igen

                    if i == 1 and j == 1 and l1 == 1:
                        icur = igen
                        p1 = F0
                        if np.float32(icur) >= rain:
                            p1 = rain
                    else:
                        if icur > 0:
                            gen = rf_amt_gen(state, gset, jc, lc, nk, p1,
                                             inext, rain, cls_y, cls_x,
                                             coffsets, ccounts, hgy, hgx,
                                             sxxi_t, sxxixy_t, scond_t, sd_t,
                                             ymax_t, wts)
                        else:
                            gen = F0

                        out[it, w] = gen
                        w += 1

                        p1 = gen
                        ip1 = np.float32(igen)
                        icur = igen

                    jp = j
                    lp = l1

        # the final day just repeats the previous value
        out[it, w] = gen
        w += 1

    return out
