"""The random stream must match the Fortran bit for bit.

Every other correctness property of the generator rests on this: the day loop
consumes one ``ran1`` per station pick and per occurrence draw, and one
``ran1`` plus one ``gasdev`` per wet day, so any drift in the stream desyncs
the whole chain.

The expected values were produced by the reference Fortran build
(``reference_run/``, gfortran ``-O2 -std=legacy``) and are the raw IEEE-754
single-precision bit patterns, so the comparison is exact rather than
approximate.
"""

import numpy as np

from pyraingen.daily.rng import gasdev, new_state, ran1, warm_up

# First values after the 500-call gasdev warm up the Fortran performs on entry.
EXPECTED_RAN1 = [
    "3F40CFFF", "3F6E00AC", "3CD7C654", "3D9E5C36", "3F16FDBF",
    "3F5B5C3C", "3D87C333", "3F6B530F", "3F018739", "3E1F24D1",
]
# ...and the deviates that follow those ten uniforms.
EXPECTED_GASDEV = [
    "BF9EF1BC", "3DB28D42", "4029EDFF", "3F2B9AB0", "3ED33147",
]


def _bits(x):
    return format(np.float32(x).view(np.uint32), "08X")


def test_ran1_matches_fortran_bit_for_bit():
    state, gset = new_state(131)
    warm_up(state, gset, 500)
    got = [_bits(ran1(state)) for _ in range(len(EXPECTED_RAN1))]
    assert got == EXPECTED_RAN1


def test_gasdev_matches_fortran_bit_for_bit():
    state, gset = new_state(131)
    warm_up(state, gset, 500)
    for _ in range(len(EXPECTED_RAN1)):
        ran1(state)
    got = [_bits(gasdev(state, gset)) for _ in range(len(EXPECTED_GASDEV))]
    assert got == EXPECTED_GASDEV


def test_ran1_stays_in_the_open_unit_interval():
    state, gset = new_state(131)
    values = np.array([ran1(state) for _ in range(20000)], dtype=np.float32)
    assert values.min() > 0.0
    assert values.max() < 1.0
    # NR caps the generator at 1 - 1.2e-7.
    assert values.max() <= np.float32(1.0) - np.float32(1.2e-7)


def test_non_negative_seeds_all_give_the_same_stream():
    """``ran1`` reseeds with ``max(-idum, 1)``, so 131 and 7 are the same run.

    This is a property of the original generator, not an oversight in the
    port, and it is why :func:`pyraingen.regionaliseddailysim` documents
    negative seeds as the way to get a different realisation set.
    """
    def first_ten(seed):
        state, gset = new_state(seed)
        warm_up(state, gset, 500)
        return [_bits(ran1(state)) for _ in range(10)]

    assert first_ten(131) == first_ten(7) == first_ten(0) == EXPECTED_RAN1


def test_negative_seeds_give_a_different_stream():
    state, gset = new_state(-99)
    warm_up(state, gset, 500)
    got = [_bits(ran1(state)) for _ in range(10)]
    assert got != EXPECTED_RAN1


def test_gasdev_caches_its_second_deviate():
    """The polar method makes two deviates per pair of uniforms."""
    state, gset = new_state(131)
    before = state[0]
    gasdev(state, gset)          # consumes uniforms, caches the partner
    mid = state[0]
    gasdev(state, gset)          # returns the cache, consumes nothing
    after = state[0]
    assert mid != before
    assert after == mid
