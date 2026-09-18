"""End-to-end regression test for the pure-Python daily generator.

The expected values come from the reference Fortran build in
``reference_run/`` -- the original ``regionalised_dailyT4.for`` with the one
change of defining ``amx``/``amn`` as the pool maximum and minimum (they were
uninitialised stack garbage before, so 1.0.2 had no well-defined output).

Run on the five bundled example stations with the historical defaults, the
port reproduces that Fortran run exactly: every one of the 153 410 generated
daily depths agrees to the 0.1 mm the Fortran wrote, and the wet/dry sequence
matches with zero flips.  The digest below is over the whole array, so any
drift in the random stream, the analogue pools or the day loop fails the test.
"""

import hashlib
from importlib import resources

import numpy as np
import pytest

from pyraingen.daily.core import run

NSIM = 10
NYEARS = 42
STARTYEAR = 1967
NDAYS = 15341

# Derived from reference_run/baseline/mmm_12345.out.
EXPECTED_SHA256 = "158ec3a0d122633c901ee8314a9ef34fdb2a381270cfedd3ab23311128da31b2"
EXPECTED_TOTAL = 521326.4
EXPECTED_WET_DAYS = 51279
EXPECTED_SIM_TOTALS = [55572.2, 51058.8, 50627.1, 56133.7, 47994.2,
                       51047.5, 47382.5, 54711.7, 56813.3, 49985.4]
EXPECTED_SIM_WET = [5223, 5054, 5164, 5457, 4957, 5016, 4882, 5246, 5277, 5003]

NEARBY = """    No Index Weight Years St_year Av annual rainfall

 target Station
     0 12345  1.000     0    -1   1100.0

 Nearby Stations
     1 61003  0.200    36  1935   1000.00
     2 61072  0.200   119  1890   1000.00
     3 61223  0.200    29  1964   1000.00
     4 66072  0.200    47  1957   1000.00
     5 66101  0.200    25  1889   1000.00
"""


@pytest.fixture(scope="module")
def example_inputs(tmp_path_factory):
    """Path to the bundled example daily data and a matching nearby file."""
    with resources.path("pyraingen.data.example.daily",
                        "rev_dr061003.txt") as f:
        data_path = str(f.parent)
    nearby = tmp_path_factory.mktemp("daily") / "nearby_station_details.out"
    nearby.write_text(NEARBY)
    return data_path, str(nearby)


@pytest.fixture(scope="module")
def generated(example_inputs):
    data_path, nearby = example_inputs
    rain, days = run(data_path=data_path, nearby_path=nearby, nyears=NYEARS,
                     startyear=STARTYEAR, nsim=NSIM, seed=131, verbose=False)
    return rain, days


def test_shape_and_calendar(generated):
    rain, days = generated
    assert rain.shape == (NDAYS, NSIM)
    assert days.shape == (NDAYS,)
    # 42 years from 1967 inclusive, with 11 leap days.
    assert NDAYS == 42 * 365 + 11
    assert np.all(np.diff(days) == 1)


def test_matches_the_fortran_reference_exactly(generated):
    rain, _ = generated
    rounded = np.round(rain.astype(np.float64), 1)
    digest = hashlib.sha256(
        np.ascontiguousarray(rounded, dtype=np.float64).tobytes()).hexdigest()
    assert digest == EXPECTED_SHA256


def test_summary_statistics(generated):
    rain, _ = generated
    rounded = np.round(rain.astype(np.float64), 1)
    assert rounded.sum() == pytest.approx(EXPECTED_TOTAL, abs=0.05)
    assert int((rounded >= 0.3).sum()) == EXPECTED_WET_DAYS
    assert rounded.sum(axis=0) == pytest.approx(EXPECTED_SIM_TOTALS, abs=0.05)
    assert [int(v) for v in (rounded >= 0.3).sum(axis=0)] == EXPECTED_SIM_WET


def test_amounts_are_physically_bounded(generated):
    rain, _ = generated
    wet = rain[rain > 0.0]
    # rf_amt_gen clamps every generated wet-day depth into [cutoff, 600].
    assert wet.min() >= np.float32(0.3)
    assert rain.max() <= np.float32(600.0)
    # Dry days are exactly zero, not a small positive residue.
    assert set(np.unique(rain[rain < np.float32(0.3)]).tolist()) == {0.0}


def test_realisations_differ(generated):
    """One RNG stream runs through all realisations, so they must not repeat."""
    rain, _ = generated
    for s in range(1, NSIM):
        assert not np.array_equal(rain[:, 0], rain[:, s])


def test_seed_is_honoured(example_inputs):
    """A negative seed selects a different stream; a non-negative one does not."""
    data_path, nearby = example_inputs
    kwargs = dict(data_path=data_path, nearby_path=nearby, nyears=3,
                  startyear=1967, nsim=1, verbose=False)
    base = run(seed=131, **kwargs)[0]
    same = run(seed=7, **kwargs)[0]
    other = run(seed=-99, **kwargs)[0]
    assert np.array_equal(base, same)
    assert not np.array_equal(base, other)
