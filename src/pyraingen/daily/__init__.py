"""Pure-Python implementation of the regionalised daily rainfall generator.

This package replaces the Fortran 77 sources that were shipped as a compiled
extension module up to pyraingen 1.0.2.  See :func:`pyraingen.daily.core.run`
for the entry point used by :func:`pyraingen.regionaliseddailysim`.
"""

from ._jit import HAVE_NUMBA

__all__ = ["HAVE_NUMBA"]
