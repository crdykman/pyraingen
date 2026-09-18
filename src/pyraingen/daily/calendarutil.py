"""Calendar walks ported from the Fortran ``day``/``daycount``/``day_neg``/
``day_pos`` routines.

Array indices follow the Fortran convention throughout this package: months
are 1..12, days 1..31 and year slots 1..nyear, with index 0 left unused.  That
keeps the ported loops a direct transliteration of the original, which matters
because several of them differ subtly from one another (notably the ``ns``
epoch offset used for leap years, which is ``nstrt-1`` in some routines and
``nstrt`` in others -- reproduced here as-is).

``nd(2)`` is 29 in the static table; callers overwrite it with
:func:`daycount` wherever a real February length is needed.
"""

import numpy as np

from ._jit import njit

NOUT = 12

# 1-based month lengths; index 0 unused.  February is 29 here, exactly as the
# Fortran ``day`` subroutine sets it.
NDAY = np.array([0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=np.int64)


@njit(cache=True)
def daycount(ns, i):
    """Length of February for calendar year ``ns + i``."""
    nd = 28
    if (ns + i) % 400 == 0:
        nd = 29
    if (ns + i) % 100 != 0 and (ns + i) % 4 == 0:
        nd = 29
    return nd


@njit(cache=True)
def month_len(nday, ns, i, j):
    """Length of month ``j`` in year slot ``i`` of a record starting at ``ns``."""
    if j == 2:
        return daycount(ns, i)
    return nday[j]


@njit(cache=True)
def day_neg(i, j, l, ip, ns, lag, nday):
    """Normalise a non-positive day-of-month by walking backwards.

    Returns ``(i, j, l, ind)`` where ``ind == 1`` means the walk ran off the
    start of the record.
    """
    if i <= ip and j == 1 and l - lag < 1:
        return i, j, l, 1

    while True:
        if l > 0:
            return i, j, l, 0
        j -= 1
        if j < 1:
            j = NOUT
            i -= 1
            if i < ip:
                return i, j, l, 1
        nd = nday[j]
        if j == 2:
            nd = daycount(ns, i)
        l += nd


@njit(cache=True)
def day_pos(i, j, l, ns, ny, nday):
    """Normalise an over-long day-of-month by walking forwards.

    Returns ``(i, j, l, ind)`` where ``ind == 1`` means the walk ran off the
    end of the record.
    """
    if i == ny and j == NOUT and l > nday[NOUT]:
        return i, j, l, 1

    while True:
        nd = nday[j]
        if j == 2:
            nd = daycount(ns, i)
        if l <= nd:
            return i, j, l, 0
        l -= nd
        j += 1
        if j > NOUT:
            i += 1
            j = 1
            if i > ny:
                return i, j, l, 1


@njit(cache=True)
def prev_day(i, j, l, ns, nday):
    """Step back one day.  Returns ``(i, j, l, ind)``; ``ind == 1`` if before
    the start of the record.

    This is the open-coded "check for prev day" block repeated in
    ``smoothprob``, ``av_sd_lon``, ``smoothavcov`` and ``rf_amt_store``.
    """
    lp = l - 1
    jp = j
    ip = i
    if lp < 1:
        if ip <= 1 and jp == 1:
            return ip, jp, lp, 1
        jp -= 1
        if jp < 1:
            ip -= 1
            jp = NOUT
        lp = nday[jp]
        if jp == 2:
            lp = daycount(ns, ip)
    return ip, jp, lp, 0


@njit(cache=True)
def iseas(j, isn):
    """Season index for month ``j`` under an ``isn``-season split."""
    if isn == 4:
        if 3 <= j <= 5:
            return 1
        if 6 <= j <= 8:
            return 2
        if 9 <= j <= 11:
            return 3
        return 4
    elif isn == 2:
        if 6 <= j <= 10:
            return 1
        return 2
    elif isn == 3:
        if 1 <= j <= 5:
            return 1
        if 6 <= j <= 9:
            return 2
        return 3
    return 1
