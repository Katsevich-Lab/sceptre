# Frozen behavior suite

These tests record what sceptre outputs on small inputs, as testthat snapshots,
and check that a change reproduces them. They are meant for changes that should
not alter results, such as the restructuring of the R internals and the move to
`SingleCellExperiment`.

The snapshots in `tests/_snaps/` record `main` after the 0.99.1 bug fixes
(#242–#245), on macOS arm64 with R 4.6.0.

This is a developer lane. It is not run by `R CMD check`, and `tools/` is
excluded from the built package (`.Rbuildignore`), so the snapshots neither
count toward the package size nor break the Bioconductor build with a
platform's last-digit differences. It runs in CI on pull requests
([`frozen-behavior.yaml`](../../.github/workflows/frozen-behavior.yaml)).

## Running

Install the version to check, then run the tests against it with the settings
CI uses:

```sh
R CMD INSTALL .
export NOT_CRAN=true TESTTHAT_MAX_FAILS=Inf LANGUAGE=en \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  VECLIB_MAXIMUM_THREADS=1 R_DATATABLE_NUM_THREADS=1
Rscript -e 'testthat::test_dir("tools/consistency/tests", package = "sceptre", load_package = "installed")'
```

The tests use whichever sceptre is installed first on the library path, so
install the version you mean to check. testthat skips snapshot tests unless
`NOT_CRAN=true`. `LANGUAGE=en` matters because errors and warnings are recorded
as text, and single threads keep results reproducible.

## Layout

- `harness.R`: the inputs, the cases and how each case is observed.
- `tests/helper-frozen.R`: runs a group of cases and snapshots each one.
- `tests/test-<group>.R`: one file per group of cases, snapshotted in
  `tests/_snaps/<group>/<case>.R` (public outputs) and
  `<case>-internal.R` (internal state).
- `tests/test-inputs.R`: the inputs, the list of cases and their test files,
  and serial against parallel.
- `tests/inputs/saved-objects/`: a saved 0.99.0 object, a file-format input.

## Cases

The cases call only exported functions. Every case starts from the same seed,
and its observation records the R random number generator state at the end.

| Case | Covers |
|---|---|
| `api` | Exports, their arguments and defaults, S4 methods, and the `sceptre_object` slots |
| `import_matrix_low` | `import_data()` from dense matrices |
| `import_matrix_vector` | `import_data()` with a `vector_id` column (guide counts summed per vector) |
| `import_dense_symmetric` | `import_data()` with a square symmetric dense matrix |
| `import_cellranger_two_directories` | `import_data_from_cellranger()` on the packaged two-directory example |
| `import_cellranger_one_directory` | The same with one directory |
| `import_cellranger_one_directory_covariates` | One directory with extra covariates; covariates only |
| `import_parse` | `import_data_from_parse()` on the packaged Parse example |
| `assign_low_maximum`, `assign_low_thresholding`, `assign_low_mixture` | `assign_grnas()` with each method on the low-MOI example data |
| `assign_high_mixture`, `assign_high_thresholding` | `assign_grnas()` on the high-MOI example data |
| `workflow_*` (7 cases) | The full workflow on synthetic data: import, parameters, assignment, QC, calibration, power and discovery. Low and high MOI; NT and complement controls; permutations and CRT; union, Bonferroni and singleton integration; each side; skew-normal or no approximation; output amounts 1 to 3 |
| `example_low_defaults` | Low-MOI example data with package defaults (maximum assignment, NT controls, permutations), and `write_outputs_to_directory()` |
| `example_high_defaults` | High-MOI example data with package defaults (mixture assignment, complement controls, CRT) |
| `pairs_constructors` | `construct_positive_control_pairs()`, `construct_cis_pairs()` and `construct_trans_pairs()` |
| `pairs_trans_caller_mutation` | Whether `construct_trans_pairs()` changes the caller's table |
| `pairs_cis_unsorted_positions` | `construct_cis_pairs()` with an unsorted position table |
| `cache_reuse` | Reuse and invalidation of cached response fits |
| `cache_control_switch` | Switching from NT to complement controls after calibration |
| `cache_control_switch_cleared` | The same with the cached fits cleared by hand |
| `edge_duplicate_pairs`, `edge_no_passing_pairs` | Rejected requests |
| `edge_intercept_only_nt` | `~ 1` with NT controls |
| `saved_object_0_99_0` | Loading a saved 0.99.0 object, and running discovery on it |

`test-inputs.R` also checks that two parallel workers give the serial results
(skipped on Windows).

**Inputs.** The synthetic data keep the guide design, covariates, pairs,
removed cells and workflow settings of the characterization suite of the
refactor (#234). The responses are negative binomial (size 2) instead of
Poisson. With Poisson counts the dispersion estimate sits near its cap of 1000
and moves with floating-point round-off, which changed p-values between
platforms. The covariate `z` is rounded to 6 digits, so the inputs are
bit-identical on every platform. The example-data cases use the data sets
shipped in `data/` and `inst/extdata/`. `test-inputs.R` snapshots the synthetic
inputs, the digests of the shipped ones and of the saved 0.99.0 object, and the
list of cases. A changed input is then reported as such rather than as changed
behavior, and a dropped case cannot silently stop being checked.

### Public outputs and internal state

Each case has two snapshots:

- **Public outputs**: what a user can see. This covers the exported getters,
  `get_result()`, the printed summary and the files of
  `write_outputs_to_directory()`. It also covers pair constructors, exports and
  their arguments, errors, warnings and the random number generator state.
- **Internal state**: `sceptre_object` slots, stored under elements named
  `internal`. These are cached response fits, resampling budgets, the QC'd cell
  and index lists, the requested pairs with their pairwise QC, and the slot
  list.

A change that should keep results must reproduce the public snapshots. A
restructuring that moves or redefines slots, such as the move to
`SingleCellExperiment`, may set `SCEPTRE_FROZEN_WAIVE_INTERNAL=true`. The cases
then read no internal slots at all, so removed or renamed slots do not break
the public snapshots, and the internal snapshots are skipped; the pull request
must say so. The cache cases are the exception: they clear the fit cache
(`response_precomputations`) as part of what they test, so a change to the
cache adapts them. Adapting
`consistency_observe()` to the new object and updating the snapshots restores
the internal checks.

## Comparison

Each snapshot is the observation written as R code with 17 significant digits,
so it reads back to within one unit in the last place and shows in a diff
which values changed. testthat
stores the snapshots, reports new and changed ones, and fails them in CI;
`consistency_differences()` in `harness.R` decides whether an observation
matches its snapshot:

- **Doubles** are compared elementwise and relatively:
  |actual − expected| ≤ tolerance · |expected|. This holds at every magnitude,
  so a p-value of 1e-25 that becomes 1e-14 is a difference. An expected value
  of exactly 0 allows |actual| ≤ 1e-12 (nothing at a tolerance of 0, which the
  serial-versus-parallel check uses). The positions of `NA`, `NaN` and `Inf`
  must match. (testthat's own `tolerance` is relative to the average size of
  the values, so it would miss the p-value above.)
- **Everything else must be identical:** integers, logicals, strings, factors
  with their levels, names, dimensions, classes, row names, row order, error
  messages, warnings, and the random number generator state after each case.

A failure lists the differing paths, for example
`observation$value$tables$discovery_result$columns$p_value`.

**Platforms.** The tolerance is `SCEPTRE_FROZEN_TOLERANCE`, 1e-8 by default.
CI runs macOS with R 4.6.0, the recording platform, at 1e-8, and it must pass.
Ubuntu runs at 1e-6 as an informational job (`continue-on-error`) until it has
passed in CI. Every decision (QC, significance, pair identity, assignments)
must be identical under both.

Known Linux difference: in `edge_intercept_only_nt`, the permutation null
statistics take few distinct values, and some equal the observed statistic.
Whether a tie counts as "at least as extreme" then depends on the last bit of
floating-point arithmetic, so the empirical p-values of `r3` and `r5` differ
by platform (0.996 on macOS, 0.94 to 0.984 on Linux). This is sceptre's
behavior, not the suite's.

CI installs `R.utils`, which `data.table::fread()` needs to read the gzipped
Cell Ranger example files.

## Updating the snapshots

When a change is meant to alter results, run the suite, then review and accept
the new snapshots:

```r
testthat::snapshot_review(path = "tools/consistency/tests")
testthat::snapshot_accept(path = "tools/consistency/tests")
```

Update the snapshots in these situations only, and explain each change in the
pull request:

- **A change is meant to alter results**, such as a bug fix.
- **The harness changed**, for example a new case or a new observed output.
- **The reference moves on**, for example to a new release.

Never update them to make an unexplained difference go away. Update on the
recording platform (macOS arm64, R 4.6.0), since other platforms can differ in
the last digits.
