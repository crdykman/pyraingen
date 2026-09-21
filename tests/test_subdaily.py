"""Tests for the subdaily disaggregation path.

These cover the defects fixed in 2.1.0, each of which survived because the
subdaily path had no assertions of any kind: its only "test" was a script that
ran the generator and printed. Several of the bugs were unreachable from the
one configuration that script used (``genSeqOption=3`` with a target absent
from the station index).

A full end-to-end run is deliberately *not* attempted. The bundled example
data does not contain a pluviograph file for every station the selection
returns -- two were missing as of 2.0.0 and six as of 2.1.0 -- so
``regionalisedsubdailysim`` cannot complete against it. The disaggregation
kernel is therefore exercised directly with synthetic pools, which also makes
these tests fast and independent of the 30 MB of bundled netCDFs.
"""

import os
from importlib import resources

import numpy as np
import pytest

from pyraingen.loadsubdailystationmeta import loadSubDailyStationMeta
from pyraingen.numberofyears import numberOfYears
from pyraingen.subdailydisaggregation import subDailyDisaggregation
from pyraingen.targetstations import targetStations

RECORDS_PER_DAY = 240
DRY_WET_CUTOFF = 0.3

# Sydney Observatory Hill -- present in index.nc, which is what makes it
# useful: the self-selection bug only fires when the target is in the index.
TARGET_INDEX = 66037
TARGET = dict(lat=-33.9410, lon=151.1730, elevation=6.00,
              distToCoast=0.27, annualRainDepth=1086.71, temperature=22.40)


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def data_paths():
    with resources.path("pyraingen.data", "index.nc") as f:
        index_path = str(f)
    with resources.path("pyraingen.data", "coefficients.dat") as f:
        coeff_path = str(f)
    return {"pathIndex": index_path, "pathCoeff": coeff_path}


@pytest.fixture(scope="module")
def station_meta(data_paths):
    return loadSubDailyStationMeta(data_paths["pathIndex"])


@pytest.fixture(scope="module")
def selection(data_paths, station_meta):
    """Run the real nearby-station selection for a target that is indexed."""
    param = {"targetIndex": TARGET_INDEX, "minYears": 10, "nYearsPool": 500}
    details, near = targetStations(param, data_paths, dict(TARGET),
                                   dict(station_meta))
    return details, near


def _numba_dict(**values):
    from numba.core import types
    from numba.typed import Dict

    out = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
    for key, value in values.items():
        out[key] = float(value)
    return out


def _fragment_pool(n_days=366, n_good=8, seed=5):
    """Build a self-consistent synthetic fragment pool.

    Consistency matters: a fragment's sub-daily profile must sum to the daily
    depth recorded for that same fragment, because the method of fragments
    rescales the chosen profile by ``target / fragmentDepth``. With unrelated
    random arrays the rescaling is meaningless and no invariant holds.
    """
    from numba.typed import List

    rng = np.random.default_rng(seed)
    fragments, state, daily_depth, n_good_days = List(), List(), List(), List()
    for _ in range(4):                       # nSeasons
        profile = rng.random((n_days, n_good, RECORDS_PER_DAY))
        profile /= profile.sum(axis=2, keepdims=True)
        depth = rng.random((n_days, n_good)) * 10 + 0.5
        fragments.append((profile * depth[:, :, None]).astype(np.float32))
        daily_depth.append(depth)
        state.append(np.ones((n_days, n_good)))
        n_good_days.append(np.full((n_days,), float(n_good)))
    return fragments, state, daily_depth, n_good_days


def _disaggregate(target_daily, year_start, year_end, seed=5):
    fragments, state, daily_depth, n_good_days = _fragment_pool(seed=seed)
    param = _numba_dict(maxNearNeighb=3, halfWinLen=5,
                        dryWetCutoff=DRY_WET_CUTOFF, absDiffTol=0.5,
                        simYearStart=year_start, simYearEnd=year_end,
                        DayStart=0, DayEnd=365)
    return subDailyDisaggregation(target_daily, param, n_good_days,
                                  fragments, state, daily_depth)


def _days_in(year_start, year_end):
    total = 0
    for year in range(year_start, year_end + 1):
        leap = (year % 400 == 0) or (year % 100 != 0 and year % 4 == 0)
        total += 366 if leap else 365
    return total


# ---------------------------------------------------------------------------
# the disaggregation kernel
# ---------------------------------------------------------------------------


def test_wet_day_totals_equal_the_daily_depth():
    """The whole point of the method: fragments are rescaled to the daily total."""
    rng = np.random.default_rng(2)
    target = rng.random(_days_in(2001, 2002)) * 5
    out = _disaggregate(target, 2001, 2002)

    wet = target >= DRY_WET_CUTOFF
    assert wet.sum() > 100, "test data should contain plenty of wet days"
    assert np.allclose(out.sum(axis=0)[wet], target[wet], rtol=1e-3)


def test_dry_days_are_exactly_zero():
    rng = np.random.default_rng(3)
    target = rng.random(_days_in(2001, 2002)) * 5
    target[::3] = 0.0
    out = _disaggregate(target, 2001, 2002)

    dry = target < DRY_WET_CUTOFF
    assert dry.any()
    assert (out[:, dry] == 0).all()


def test_output_is_physically_valid():
    rng = np.random.default_rng(4)
    target = rng.random(_days_in(2001, 2001)) * 5
    out = _disaggregate(target, 2001, 2001)

    assert out.shape == (RECORDS_PER_DAY, target.size)
    assert np.isfinite(out).all()
    assert (out >= 0).all()


@pytest.mark.parametrize("year_start,year_end,label", [
    (2001, 2004, "final year is a leap year"),
    (2001, 2003, "final year is not a leap year"),
])
def test_final_days_are_disaggregated_correctly(year_start, year_end, label):
    """Regression test for the end-of-sequence branch.

    It used to test ``loopDay == nDaysCurrYear-1``, but the day loop always
    runs over 366 slots and skips 29 February in non-leap years, so the last
    day sits at index 365 either way. In a non-leap final year the branch
    therefore fired a day early *and* again at 365, and because it never
    shifted the day window both of those days were disaggregated against the
    previous day's depth. Checking the total on the last few days catches it.
    """
    rng = np.random.default_rng(7)
    target = rng.random(_days_in(year_start, year_end)) * 5 + 1.0   # all wet
    out = _disaggregate(target, year_start, year_end)

    assert out.shape[1] == target.size, label
    totals = out.sum(axis=0)
    assert np.allclose(totals[-3:], target[-3:], rtol=1e-3), (
        f"{label}: last three days should honour their own daily depth, "
        f"got {totals[-3:]} for {target[-3:]}"
    )


def test_kernel_releases_the_gil():
    """regionalisedsubdailysim dispatches this through joblib's thread pool.

    A numba nopython function holds the GIL unless told otherwise, so without
    nogil the "parallel" simulations run strictly one at a time.
    """
    assert subDailyDisaggregation.targetoptions.get("nogil") is True


def test_kernel_does_not_mutate_its_shared_inputs():
    """Threads share the fragment pools, so the kernel must only read them."""
    fragments, state, daily_depth, n_good_days = _fragment_pool()
    before = [np.array(fragments[s], copy=True) for s in range(4)]
    depth_before = [np.array(daily_depth[s], copy=True) for s in range(4)]

    param = _numba_dict(maxNearNeighb=3, halfWinLen=5,
                        dryWetCutoff=DRY_WET_CUTOFF, absDiffTol=0.5,
                        simYearStart=2001, simYearEnd=2001, DayStart=0,
                        DayEnd=365)
    target = np.random.default_rng(8).random(365) * 5
    subDailyDisaggregation(target, param, n_good_days, fragments, state,
                           daily_depth)

    for s in range(4):
        assert np.array_equal(fragments[s], before[s])
        assert np.array_equal(daily_depth[s], depth_before[s])


# ---------------------------------------------------------------------------
# nearby-station selection
# ---------------------------------------------------------------------------


def test_target_is_not_its_own_nearest_neighbour(selection, station_meta):
    """``np.where`` returns a tuple, so the self-exclusion guard never fired.

    The target then scored a perfect similarity against itself and ranked
    first in its own nearby-station list.
    """
    details, near = selection
    target_pos = int(np.flatnonzero(
        np.asarray(details["stnIndex"]) == TARGET_INDEX)[0])

    for season in range(near.shape[0]):
        chosen = [int(i) for i in near[season, :] if i >= 0]
        assert target_pos not in chosen, (
            f"season {season} selected the target station as a neighbour"
        )


def test_unused_slots_use_a_negative_sentinel(selection):
    """0 is a valid station index, so it cannot mean "no more stations"."""
    _details, near = selection
    unused = near[near < 0]
    assert unused.size > 0, "expected some unfilled slots"
    assert (unused == -1).all()


def test_selected_indices_are_in_range(selection):
    details, near = selection
    n_stations = len(details["stnIndex"])
    chosen = near[near >= 0]
    assert chosen.size > 0
    assert chosen.max() < n_stations
    assert (chosen == chosen.astype(int)).all()


def test_station_index_zero_is_not_treated_as_a_terminator(data_paths,
                                                           station_meta):
    """Consumer side of the sentinel fix.

    ``dailySequences``, ``getFragments`` and ``numberOfYears`` all break out of
    the station loop on the sentinel. While that sentinel was 0, a season whose
    most similar station happened to be first in the metadata silently got an
    empty pool. Here station 0 is a genuine selection and must be counted.
    """
    meta = dict(station_meta)
    # Point slot 0 at a station that has a bundled pluviograph file.
    with resources.path("pyraingen.data.example.subdaily", "daily.nc") as f:
        bundled_dir = str(f.parent)
    available = sorted(int(n[3:9]) for n in os.listdir(bundled_dir)
                       if n.startswith("plv"))
    assert available, "expected bundled pluviograph files"

    index = np.asarray(meta["stnIndex"]).copy()
    index[0] = available[0]
    meta["stnIndex"] = index

    near = np.full((4, 3), -1.0)
    near[:, 0] = 0                       # station 0 is the only selection

    paths = dict(data_paths)
    paths["pathSubDaily"] = bundled_dir
    years = numberOfYears(4, meta, near, paths)

    assert (np.asarray(years) > 0).all(), (
        "station index 0 was skipped, leaving an empty pool"
    )


# ---------------------------------------------------------------------------
# regionalisedsubdailysim entry point
# ---------------------------------------------------------------------------


def test_unknown_target_is_rejected_clearly(tmp_path):
    from pyraingen.regionalisedsubdailysim import regionalisedsubdailysim

    with pytest.raises(ValueError, match="Target station not found"):
        regionalisedsubdailysim(
            str(tmp_path / "missing.nc"), str(tmp_path) + os.sep, 999999,
            pathIndex=None, pathCoeff=None, pathReference=" ",
            fnameSubDaily=str(tmp_path / "out.nc"),
            minYears=10, nYearsPool=50, dryWetCutoff=DRY_WET_CUTOFF,
            halfWinLen=15, maxNearNeighb=10, nSims=1,
            genSeqOption=0, nYearsRef=5, absDiffTol=0.1,
        )


def test_known_target_gets_past_the_index_lookup(tmp_path):
    """genSeqOption 0, 1 and 4 used to die on ``ndarray.index()``.

    Those modes could never run: the only configurations exercised anywhere
    were 2 and 3, which take the other branch. The call below still fails --
    there is no input data -- but it must get past the lookup, so an
    AttributeError here is the regression.
    """
    from pyraingen.regionalisedsubdailysim import regionalisedsubdailysim

    with pytest.raises(Exception) as excinfo:
        regionalisedsubdailysim(
            str(tmp_path / "missing.nc"), str(tmp_path) + os.sep, TARGET_INDEX,
            pathIndex=None, pathCoeff=None, pathReference=" ",
            fnameSubDaily=str(tmp_path / "out.nc"),
            minYears=10, nYearsPool=50, dryWetCutoff=DRY_WET_CUTOFF,
            halfWinLen=15, maxNearNeighb=10, nSims=1,
            genSeqOption=0, nYearsRef=5, absDiffTol=0.1,
        )

    assert not isinstance(excinfo.value, AttributeError), (
        f"the ndarray .index() regression is back: {excinfo.value}"
    )
    assert "index" not in str(excinfo.value) or "No such file" in str(excinfo.value)
