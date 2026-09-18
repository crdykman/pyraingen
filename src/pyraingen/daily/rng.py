"""Bit-exact ports of the Numerical Recipes generators used by the Fortran.

``ran1`` is the NR minimal-standard generator with a Bays-Durham shuffle and
``gasdev`` is the polar Box-Muller transform with the cached second deviate.
Both are pure integer / float32 arithmetic, so the stream is reproducible
across platforms and matches the Fortran exactly.

The generator state lives in caller-owned arrays so it can be threaded
through ``@njit`` kernels:

* ``state`` -- ``int64[35]``: ``idum``, ``iy``, ``iset``, then the 32-entry
  shuffle table ``iv``.
* ``gset`` -- ``float32[1]``: the cached second normal deviate.

.. note::
   ``ran1`` reseeds with ``idum = max(-idum, 1)``.  The Fortran calls it with
   ``iseed = 131``, so the stream actually starts from ``idum = 1``: *every*
   non-negative seed gives the identical stream.  Pass a negative seed to
   select a different stream.
"""

import math

import numpy as np

from ._jit import njit

NTAB = 32
IA = 16807
IM = 2147483647
IQ = 127773
IR = 2836
NDIV = 1 + (IM - 1) // NTAB

# Both constants are formed in single precision exactly as the Fortran does.
AM = np.float32(1.0) / np.float32(IM)
RNMX = np.float32(1.0) - np.float32(1.2e-7)

STATE_SIZE = 3 + NTAB

_F0 = np.float32(0.0)
_F1 = np.float32(1.0)
_F2 = np.float32(2.0)
_FM2 = np.float32(-2.0)


def new_state(seed=131):
    """Allocate a fresh generator state for ``seed``.

    Parameters
    ----------
    seed : int
        Seed passed to ``ran1``.  Matching the Fortran, any value ``>= 0``
        yields the same stream; use a negative value for a distinct stream.

    Returns
    -------
    tuple of ndarray
        ``(state, gset)`` as described in the module docstring.
    """
    state = np.zeros(STATE_SIZE, dtype=np.int64)
    state[0] = int(seed)
    gset = np.zeros(1, dtype=np.float32)
    return state, gset


@njit(cache=True)
def ran1(state):
    """Return the next uniform deviate on (0, 1) and advance ``state``."""
    idum = state[0]
    iy = state[1]

    if idum <= 0 or iy == 0:
        idum = -idum
        if idum < 1:
            idum = 1
        for j in range(NTAB + 8, 0, -1):
            k = idum // IQ
            idum = IA * (idum - k * IQ) - IR * k
            if idum < 0:
                idum += IM
            if j <= NTAB:
                state[2 + j] = idum
        iy = state[3]

    k = idum // IQ
    idum = IA * (idum - k * IQ) - IR * k
    if idum < 0:
        idum += IM
    j = 1 + iy // NDIV
    iy = state[2 + j]
    state[2 + j] = idum
    state[0] = idum
    state[1] = iy

    value = AM * np.float32(iy)
    if value > RNMX:
        value = RNMX
    return value


@njit(cache=True)
def gasdev(state, gset):
    """Return the next standard normal deviate and advance ``state``."""
    if state[2] < 1:
        while True:
            v1 = _F2 * ran1(state) - _F1
            v2 = _F2 * ran1(state) - _F1
            rsq = v1 * v1 + v2 * v2
            if rsq < _F1 and rsq != _F0:
                break
        # math.log/math.sqrt evaluate in float64 and round once, which gives
        # the same answer compiled or interpreted.  numpy's float32 scalar
        # log differs from numba's by an ulp, which would desync the two
        # paths; this spelling matches the Fortran in both.
        fac = np.float32(math.sqrt(_FM2 * np.float32(math.log(rsq)) / rsq))
        gset[0] = v1 * fac
        state[2] = 1
        return v2 * fac

    state[2] = 0
    return gset[0]


@njit(cache=True)
def warm_up(state, gset, ncalls):
    """Burn ``ncalls`` normal deviates, as the Fortran does on entry."""
    for _ in range(ncalls):
        gasdev(state, gset)
