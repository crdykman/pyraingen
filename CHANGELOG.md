# Changelog

## v2.1.0 (18/09/2026)

Correctness and performance work across the subdaily and IFD code, from a
review of the whole package. **Generated results change** for the subdaily
path and for daily runs whose target is an existing gauge; the daily
generator's reproduction of the reference Fortran is unaffected and its test
suite still passes unchanged.

Supersedes 2.0.1, which was tagged but never released to PyPI.

### Fixed

- **`ifdcond` could not be imported on numpy >= 2.** `np.RankWarning` was
  removed from the top-level numpy namespace in 2.0, so `import` raised
  `AttributeError` and one of the three headline functions was unusable in a
  fresh 2.0.0 install. This was a regression introduced in 2.0.0 by relaxing
  the `numpy==1.23.5` pin.
- **`genSeqOption` 0, 1 and 4 always failed.** `regionalisedsubdailysim` called
  the list method `.index()` on an ndarray. Only `genSeqOption` 2 and 3 took
  the other branch, and both the test and the documented example use 3, so
  three of the five modes had never run.
- **Nearby-station search discarded up to 45% of candidates.** The exclusion
  test in `getnearbystations` omitted the `abs()` calls its own comment
  specifies, so it excluded every station both west *and* south of the target.
  For one real gauge it removed 305 of 1694 stations instead of 1. Only
  applies when the target index is itself in `stn_record.csv`.
- **The subdaily target station selected itself as its own nearest
  neighbour.** `np.where()` returns a tuple, so the `loopStation != idxTarget`
  guard in `targetstations` was never false and the self-exclusion the code
  intends never happened. This was also the root cause of the numpy 2
  `TypeError` patched in 2.0.0, whose workaround is now removed.
- **`0` was used as a "no more stations" sentinel** in `nearStationIdx`, but
  `0` is a valid index. If the most similar station was the first in the
  metadata, that season silently got an empty pool; any station after index 0
  was dropped. The sentinel is now `-1`, across all five consumers.
- **The subdaily end-of-sequence branch fired twice** when the final simulated
  year was not a leap year, disaggregating the last two days as if each were
  the end of the record. It now also shifts the day window, so the final day
  uses its own depth rather than the previous day's.
- **`computeIFD` and `aggregateRainfall` ignored `nRecordsPerDay` and
  `minsPerSample`**, hardcoding 240 and 6, so any non-six-minute input was
  silently mis-aggregated despite `ifdcond` exposing both as parameters.
  Aggregation now also refuses a duration that does not divide evenly rather
  than discarding the remainder.
- **A single-duration IFD run crashed twice over**: `plt.subplots` returns a
  bare Axes for one row, and `np.loadtxt` returns a 1-D array for a
  one-column target IFD. Both are fixed, and `ifdcondobjfun` no longer breaks
  at `nSims == 1`.
- `scipy.optimize.fmin` convergence is now checked and warned about instead of
  silently accepted, and no longer prints to stdout.
- Seven netCDF reader handles were never closed, three of them inside loops
  over seasons and stations.

### Changed

- `ifdcond` takes a `seed` argument, making its subsample of simulations
  reproducible. Default `None` preserves current behaviour.
- `ifdcond` now raises a clear error when `TargetIFDdurations` differs from
  `TargetIFDdurationsEst`. The docstring described a proper subset as
  supported, but `ifdcondobjfun` indexes the simulated IFD with a counter over
  the target's columns, so a subset raised `IndexError`. The limitation is now
  documented and enforced rather than hit by surprise.
- `numberOfYears` documents that it counts calendar years *spanned*, not years
  of usable record; gaps and `stnIPer` are ignored. Behaviour unchanged.
- Removed `subdailysimloop0.py` and `subdailysimloop1.py`, referenced only
  from commented-out imports.

### Performance

| change | effect |
| ------ | ------ |
| `getnearbystations.station()` vectorised | **208 s -> 0.033 s**, output identical (max difference 0.000e+00) |
| IFD annual-maximum loop slices rows instead of copying | ~**400x** on that loop; avoids ~99 GB of copies at `nSims=10` |
| `subDailyDisaggregation` marked `nogil=True` | joblib's thread pool actually parallelises: measured 0.92x -> **2.89x** on 4 threads |
| `fragments` stored as float32 | peak allocation ~2.8 GB -> ~1.4 GB; error 2.4e-5 mm, 4096x finer than the 0.1 mm data |
| redundant boolean copy dropped in IFD step 2 | ~12% faster, 37 MB less per recursion and duration |

Before this release `station()` dominated a daily run: the generator itself
takes 1.7 s while the nearby-station search took over three minutes.

## v2.0.1 (18/09/2026)

Packaging and documentation housekeeping. No functional change: the wheel's
payload is byte-for-byte the same set of files as 2.0.0.

- **Read the Docs references removed**, as the site is no longer available for
  this package. Deletes `.readthedocs.yml`, drops the
  `pyraingen.readthedocs.io` link from the README in favour of pointing at
  `docs/example.ipynb`, and removes the `sphinx-rtd-theme` dependency. The
  Sphinx config now uses `alabaster`, which ships with Sphinx, so `docs/` still
  builds locally.
- **Build backend switched from `poetry-core` to `setuptools`**, with metadata
  moved to the standard PEP 621 `[project]` table. `python -m build` remains
  the way to build. `poetry.lock` is removed.
  - Package data is now declared explicitly via
    `[tool.setuptools.package-data]`. setuptools ships only `.py` files by
    default, and the generator loads its station record, regression
    coefficients and example data through `importlib.resources`, so the
    `.csv`, `.dat`, `.nc` and `.txt` files have to be named. The built wheel
    was diffed against the published 2.0.0 wheel to confirm the file set is
    unchanged.
- The `license` metadata field is **no longer set**. poetry recorded the
  literal string `"None"`, which is not a valid SPDX expression and is
  rejected by `setuptools>=77`. The terms remain stated in `README.md`; set an
  SPDX expression in `pyproject.toml` if the metadata should carry them.

## v2.0.0 (18/09/2026)

The regionalised daily generator is now pure Python. The Fortran 77
implementation and its compiled extension module have been removed.

### Removed

- `src/pyraingen/fortran_daily/` and the two `cp38-win_amd64.pyd` binaries.
  The Fortran sources are kept for reference in `legacy/fortran/` but are no
  longer part of the installed package.
- `get_fortran_data.copy_fortran_data()` and `get_for_path()`. Nothing is
  copied into the working directory any more.
- The `data_r.dat` round trip. Parameters are passed as arguments, and the
  generated series goes straight into the netCDF instead of via a
  `mmm_<idx>.out` text file. `convertdailync()` is still available for anyone
  with existing text output.

### Changed

- **The package is no longer restricted to Python 3.8.5 on Windows.** The
  wheel is `py3-none-any`; dependency pins are now lower bounds rather than
  exact versions.
- `regionaliseddailysim()` keeps its signature and adds a `seed` keyword
  (default 131, the value the Fortran hard-coded).
- `output_stats` is accepted but ignored. It never had any effect: the Fortran
  read the argument and then wrote to a hard-coded `stat.out`. No diagnostics
  file is written now.
- Bias correction uses `groupby("day.year")` instead of `resample(day="A")`,
  which pandas renamed in 2.2. The computed scaling is unchanged.
- `numba` remains a dependency and is what makes the generator fast, but the
  daily generator now degrades gracefully to plain Python without it, and the
  two paths produce bit-identical output.

### Fixed

- `targetstations.py` raised `TypeError: only 0-dimensional arrays can be
  converted to Python scalars` on numpy >= 2, which broke the *subdaily*
  nearby-station selection. When the target site is already in the station
  list, `idxTarget` comes from `np.where()` and is a tuple, so each attribute
  delta is a one-element array; numpy 1.x converted that silently for
  `math.exp()` and numpy 2 does not. The value is now converted explicitly,
  which reproduces the previous behaviour in both cases. This was latent
  before only because the package pinned `numpy==1.23.5`.

### Results

Occurrences are unchanged from 1.0.2. The wet/dry sequence for a given seed is
bit-for-bit what the Fortran produced, because the `ran1`/`gasdev` streams are
exact integer/float32 ports and each wet day consumes the same number of
draws.

Against the reference Fortran build, over **100 realisations x 42 years on the
five bundled example stations -- 1 534 100 daily values -- the two outputs are
identical**, with no cell differing and no wet/dry flips.

**Wet-day amounts do change**, because of the `amx`/`amn` fix below. Measured
over 10 realisations x 42 years on the five bundled example stations:

| statistic              | change |
| ---------------------- | -----: |
| wet-day frequency      | 0.000% |
| mean wet-day depth     | +1.34% |
| SD of wet-day depth    | +3.16% |
| mean annual total      | +1.34% |
| mean annual maximum    | +2.21% |
| lag-1 autocorrelation  | -4.20% |
| wet/dry sequence       | identical |

- **`amx`/`amn` fix.** In `hfracx` and `hfracy` the +/-0.5*IQR window around
  each pool point was clipped against `amx` and `amn`, which were never
  assigned anywhere in the Fortran. They were uninitialised stack values, so
  bandwidths near the pool extremes depended on the compiler and on call
  history, and 1.0.2 had no single well-defined output. They are now the pool
  maximum and minimum, which is the evident intent. The percentages above are
  against a gfortran build of the original; the shipped `.pyd` was a different
  build and may have differed again.

- Two further oddities in the Fortran are **reproduced deliberately**, because
  changing them would alter results for no clear benefit. Both are documented
  in the source:
  - `simulate` guards the `(2*pi)^(-nlon/2) * det^(-1/2)` normalisation of the
    two multivariate-normal densities with `det`/`ddet`, which are declared but
    never assigned, so the normalisation never applies.
  - `hfracx` recomputes its bandwidth without the 0.8 factor that `psimmain`
    applies, so the kernel weights and the per-point bandwidths use different
    scalings.

- `seed` values `>= 0` all give the same stream. `ran1` reseeds with
  `max(-seed, 1)`, so the historical default of 131 actually started the
  generator from 1. Pass a **negative** seed for a different realisation set.

### Performance

On the bundled example stations, 42 years, single core:

| build | 10 realisations | 100 realisations |
| ----- | --------------: | ---------------: |
| Fortran (gfortran -O2) | 13.6 s | 37.1 s |
| Python + numba, warm cache | 1.7 s | 2.7 s |
| Python + numba, first call (includes JIT) | 13.3 s | 15.7 s |

The gap widens with the number of realisations because the one-off JIT and
pool-statistics costs amortise.

The gain comes from hoisting `rank_h`, `hfracx`, `hfracy` and the conditional
covariance block out of the day loop: they depend only on the analogue pool,
not on the day, and the Fortran recomputed them for all 79 157 wet days.
Without numba the generator is correct but far slower, and is not recommended
for production runs.

## v1.0.2 (18/12/2022)

- Fix bug in get get_example_data.py

## v1.0.1 (17/12/2022)

- set dependencies in pyproject.toml as static 

## v1.0.0 (16/12/2022)

- Major update to regionalised daily sim nearby site selection.
- Minor bug fixes
- Updated example
- Updated example data

## v0.2.3 (01/12/2022)

- ????

## v0.2.2 (30/11/2022)

- Includes example datasets.

## v0.1.0 (21/11/2022)

- First release of `pyraingen`!
