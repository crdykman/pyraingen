"""Optional numba acceleration.

The kernels in this package are written in the numba-compatible subset of
Python (typed arrays, explicit loops, no Python objects) so that the same
source runs either compiled or interpreted.  If numba is not installed,
:func:`njit` degrades to a no-op decorator and everything still works, just
slower.

Set ``NUMBA_DISABLE_JIT=1`` in the environment to exercise the interpreted
path while numba is installed; numba honours that itself.
"""

try:  # pragma: no cover - trivial import guard
    from numba import njit as _numba_njit

    HAVE_NUMBA = True
except ImportError:  # pragma: no cover - exercised only without numba
    HAVE_NUMBA = False


if HAVE_NUMBA:
    njit = _numba_njit
else:  # pragma: no cover - exercised only without numba

    def njit(*args, **kwargs):
        """No-op stand-in for :func:`numba.njit`.

        Accepts both the bare ``@njit`` and the parameterised
        ``@njit(cache=True)`` spellings.
        """
        if len(args) == 1 and callable(args[0]) and not kwargs:
            return args[0]

        def _decorate(func):
            return func

        return _decorate


__all__ = ["njit", "HAVE_NUMBA"]
