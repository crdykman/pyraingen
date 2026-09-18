"""Kernel-level tests, including the no-numba fallback.

The kernels are written in the numba-compatible subset so the same source runs
compiled or interpreted.  Rather than running the whole generator twice (the
interpreted path is orders of magnitude slower), these tests call each
kernel's underlying Python function via ``.py_func`` and check it agrees with
the compiled one.

Agreement is exact, and staying that way needs care.  ``np.exp`` and
``np.log`` on float32 *scalars* round differently in numba's compiled code
than in numpy's, which desynced the two paths and left ~1% of wet-day depths
differing in their last bit.  The kernels therefore call ``math.exp`` /
``math.log`` / ``math.sqrt``, which evaluate in float64 and round once, giving
the same answer either way -- and the same answer as the Fortran.
"""

import numpy as np
import pytest

from pyraingen.daily import preprocess
from pyraingen.daily.calendarutil import NDAY, day_neg, day_pos, daycount, iseas
from pyraingen.daily.poolstats import hfrac, rank_h
from pyraingen.daily.rng import gasdev, new_state, ran1, warm_up


def interpreted(fn):
    """The plain-Python version of a kernel, whether or not numba is present."""
    return getattr(fn, "py_func", fn)


# --------------------------------------------------------------------------
# calendar
# --------------------------------------------------------------------------


@pytest.mark.parametrize("year,expected", [
    (1900, 28), (1996, 29), (1999, 28), (2000, 29), (2004, 29), (2100, 28),
])
def test_daycount_leap_rule(year, expected):
    assert daycount(year, 0) == expected


def test_day_pos_rolls_february_by_leap_year():
    # 1968 is a leap year: slot 2 of a record starting 1967 keeps 29 February.
    i, j, l, ind = day_pos(2, 2, 29, 1966, 50, NDAY)
    assert (i, j, l, ind) == (2, 2, 29, 0)
    # 1969 is not: the same day rolls into 1 March.
    i, j, l, ind = day_pos(3, 2, 29, 1966, 50, NDAY)
    assert (i, j, l, ind) == (3, 3, 1, 0)


def test_day_pos_flags_running_off_the_end():
    _i, _j, _l, ind = day_pos(50, 12, 32, 1966, 50, NDAY)
    assert ind == 1


def test_day_neg_walks_back_across_a_year_boundary():
    i, j, l, ind = day_neg(5, 1, -2, 1, 1966, 2, NDAY)
    assert ind == 0
    assert (i, j) == (4, 12)
    assert l == 29          # 31 December minus two days


def test_day_neg_flags_running_off_the_start():
    _i, _j, _l, ind = day_neg(1, 1, 1, 1, 1966, 2, NDAY)
    assert ind == 1


def test_iseas_covers_every_month():
    assert [iseas(m, 4) for m in range(1, 13)] == [4, 4, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4]


# --------------------------------------------------------------------------
# solve
# --------------------------------------------------------------------------


def _matmul3(a, b):
    """3x3 product without BLAS, so the test does not depend on the backend."""
    return (a[:, :, None] * b[None, :, :]).sum(axis=1)


def _det3(m):
    return float(
        m[0, 0] * (m[1, 1] * m[2, 2] - m[1, 2] * m[2, 1])
        - m[0, 1] * (m[1, 0] * m[2, 2] - m[1, 2] * m[2, 0])
        + m[0, 2] * (m[1, 0] * m[2, 1] - m[1, 1] * m[2, 0])
    )


def test_solve_inverts_a_symmetric_matrix():
    a = np.array([[4.0, 1.0, 0.5], [1.0, 3.0, 1.0], [0.5, 1.0, 2.0]],
                 dtype=np.float32)
    work = a.copy()
    det = preprocess.solve(work, 3)
    identity = _matmul3(work.astype(np.float64), a.astype(np.float64))
    assert np.allclose(identity, np.eye(3), atol=1e-4)
    assert det == pytest.approx(abs(_det3(a.astype(np.float64))), rel=1e-4)


def test_solve_scalar_case_leaves_a_zero_variance_alone():
    """``solve`` with nv == 1 must not invert zero -- downstream relies on it."""
    work = np.zeros((1, 1), dtype=np.float32)
    det = preprocess.solve(work, 1)
    assert det == 0.0
    assert work[0, 0] == 0.0

    work = np.array([[4.0]], dtype=np.float32)
    det = preprocess.solve(work, 1)
    assert det == pytest.approx(4.0)
    assert work[0, 0] == pytest.approx(0.25)


def test_solve_matches_the_interpreted_path():
    a = np.array([[6.0, 2.0, 1.0], [2.0, 5.0, 2.0], [1.0, 2.0, 4.0]],
                 dtype=np.float32)
    fast, slow = a.copy(), a.copy()
    d1 = preprocess.solve(fast, 3)
    d2 = interpreted(preprocess.solve)(slow, 3)
    assert np.allclose(fast, slow, rtol=1e-6, atol=1e-7)
    assert d1 == pytest.approx(d2, rel=1e-6)


# --------------------------------------------------------------------------
# pool statistics
# --------------------------------------------------------------------------


def test_rank_h_is_the_quartile_gap():
    v = np.arange(1.0, 101.0, dtype=np.float32)
    # 1-based indices int(100*0.25+0.05)=25 and int(100*0.75+0.05)=75 into a
    # descending sort, i.e. the values 76 and 26.
    assert rank_h(v, 100) == pytest.approx(50.0)


def test_rank_h_guards_tiny_pools():
    """int(n*0.25+0.05) is 0 for small n, which would index out of bounds."""
    for n in (1, 2, 3, 5):
        v = np.arange(1.0, n + 1.0, dtype=np.float32)
        assert np.isfinite(rank_h(v, n))


def test_hfrac_stays_inside_its_clamps():
    rng = np.random.default_rng(1)
    v = np.sort(rng.gamma(2.0, 5.0, 200).astype(np.float32))[::-1].copy()
    out = np.zeros(200, dtype=np.float32)
    hfrac(v, 200, out)
    href = np.float32((4.0 / 3.0) ** 0.2) * np.float32(200.0 ** -0.2)
    assert out.min() >= href / np.float32(3.5) * np.float32(0.999)
    assert out.max() <= href * np.float32(3.5) * np.float32(1.001)


def test_hfrac_matches_the_interpreted_path():
    rng = np.random.default_rng(7)
    v = np.sort(rng.gamma(2.0, 5.0, 150).astype(np.float32))[::-1].copy()
    fast = np.zeros(150, dtype=np.float32)
    slow = np.zeros(150, dtype=np.float32)
    hfrac(v, 150, fast)
    interpreted(hfrac)(v, 150, slow)
    assert np.allclose(fast, slow, rtol=1e-6, atol=1e-7)


# --------------------------------------------------------------------------
# rank1: the unstable sort has to stay unstable
# --------------------------------------------------------------------------


def test_rank1_sorts_descending_and_carries_its_companions():
    xc = np.array([1.0, 3.0, 2.0, 3.0], dtype=np.float32)
    xp = np.array([10.0, 30.0, 20.0, 31.0], dtype=np.float32)
    ixn = np.array([0, 1, 0, 1], dtype=np.int64)
    preprocess.rank1(xc, xp, ixn, 4)
    assert list(xc) == [3.0, 3.0, 2.0, 1.0]
    # companions travel with their value
    assert list(xp) == [30.0, 31.0, 20.0, 10.0]
    assert list(ixn) == [1, 1, 0, 0]


def test_rank1_is_not_a_stable_sort():
    """Documents the tie order the generated series actually depends on.

    A stable descending sort of [1, 2, 1, 2] would give companions
    [b, d, a, c]; the Fortran's swap-based sort gives [b, d, c, a].  Replacing
    this with ``np.argsort(kind="stable")`` changes which analogue is drawn and
    so changes the output.
    """
    xc = np.array([1.0, 2.0, 1.0, 2.0], dtype=np.float32)
    xp = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float32)   # a, b, c, d
    ixn = np.zeros(4, dtype=np.int64)
    preprocess.rank1(xc, xp, ixn, 4)
    assert list(xc) == [2.0, 2.0, 1.0, 1.0]
    assert list(xp) == [2.0, 4.0, 3.0, 1.0]                 # b, d, c, a


# --------------------------------------------------------------------------
# rng under both paths
# --------------------------------------------------------------------------


def test_rng_is_identical_interpreted_and_compiled():
    a_state, a_gset = new_state(131)
    b_state, b_gset = new_state(131)
    warm_up(a_state, a_gset, 500)
    interpreted(warm_up)(b_state, b_gset, 500)
    assert np.array_equal(a_state, b_state)

    for _ in range(200):
        assert ran1(a_state) == interpreted(ran1)(b_state)
    for _ in range(200):
        assert gasdev(a_state, a_gset) == interpreted(gasdev)(b_state, b_gset)


def test_basic_matches_numpy_moments():
    v = np.array([1.0, 2.0, 4.0, 8.0, 16.0], dtype=np.float32)
    ave, sd = preprocess.basic(v, 5)
    assert ave == pytest.approx(6.2, rel=1e-5)
    assert sd == pytest.approx(np.std(v.astype(np.float64), ddof=1), rel=1e-4)


def test_basic_handles_degenerate_pools():
    v = np.array([3.0], dtype=np.float32)
    ave, sd = preprocess.basic(v, 1)
    assert (ave, sd) == (0.0, 1.0)
