# Frozen behavior suite

These tests record what a reference commit of sceptre outputs on small inputs,
and check that a candidate reproduces it. They are meant for changes that should
not alter results, such as the restructuring of the R internals and the move to
`SingleCellExperiment`. The fixtures are also the R reference that pysceptre is
meant to match ([fixture format](#fixture-format)).

The current fixtures record `main` at
`21f9ea098b69c444a884a49bc254de711968e3b5` (sceptre 0.99.0), on macOS arm64
with R 4.6.0. The reference commit is a parameter:
[regenerating](#regenerating-the-fixtures) records any commit.

This is a developer lane. It is not run by `R CMD check`, and `tools/` is
excluded from the built package (`.Rbuildignore`). The reasons:

- The fixtures are outputs recorded on one platform. A numerical difference on
  another platform must not break the Bioconductor build.
- The fixtures take 1.7 MB, which would count toward the package size.

The lane runs in CI on pull requests
([`frozen-behavior.yaml`](../../.github/workflows/frozen-behavior.yaml)).

## Running

You need a C++ toolchain, sceptre's dependencies, and `testthat`, `withr` and
`jsonlite`. R 4.5 or later provides `tools::sha256sum()`; on older R, install
`digest`. From the repository root:

```sh
# Check the working tree: tracked files, with their uncommitted changes
Rscript tools/consistency/consistency.R run

# Check a committed ref, or a library in which sceptre is already installed
Rscript tools/consistency/consistency.R run --ref my-branch
Rscript tools/consistency/consistency.R run --lib /path/to/library

# Run only some test files (a regular expression), or choose the tier
Rscript tools/consistency/consistency.R run --filter public
Rscript tools/consistency/consistency.R run --tier strict
```

`run` installs the candidate into a temporary library, then runs
[`testthat/`](testthat) against it in a fresh R process, with messages in
English. Untracked files are not installed: `git add` a new source file first.
Everything `run` writes goes to `tempdir()`. A full run takes about 35 seconds
on an M-series Mac, of which about 20 seconds is the installation.

To compare two commits without using the stored fixtures, record the reference
into a temporary directory and check the candidate against it:

```sh
Rscript tools/consistency/consistency.R compare --reference main --ref my-branch
```

## Layout

| Path | Contents |
|---|---|
| `harness.R` | Inputs, cases and how each case is observed. Used both to record and to test. |
| `fixture-io.R` | Digests of files and inputs. |
| `export.R` | The language-neutral copy: TSV tables and Matrix Market inputs. |
| `consistency.R` | The command: `run`, `regenerate`, `compare` and `record`. |
| `testthat/` | Public outputs and internal state of each case (separate tests), a two-worker run, fixture integrity, the TSV copy, and checks that the comparison detects small changes. |
| `fixtures/manifest.json` | Provenance, comparison rules, input digests, known-bug flags, and a description of every case and table. |
| `fixtures/cases/<case>.rds` | The full observation of each case. |
| `fixtures/tables/<case>/<table>.tsv` | The case's tables, in a language-neutral form. |
| `fixtures/inputs/` | Language-neutral copies of the inputs. |
| `fixtures/saved-objects/` | A saved object written by sceptre 0.99.0. It is an input and is never regenerated. |

## Cases

The cases call only exported functions. Every case starts from the same seed,
and its observation records the R random number generator state at the end.

| Case | Covers |
|---|---|
| `api` | Exports, their arguments and defaults, S4 methods, and the `sceptre_object` slots |
| `import_matrix_low` | `import_data()` from dense matrices |
| `import_matrix_vector` | `import_data()` with a `vector_id` column (guide counts summed per vector) |
| `import_dense_symmetric` | `import_data()` with a square symmetric dense matrix (known bug) |
| `import_cellranger_two_directories` | `import_data_from_cellranger()` on the packaged two-directory example |
| `import_cellranger_one_directory` | The same with one directory |
| `import_cellranger_one_directory_covariates` | One directory with extra covariates; covariates only (known bug) |
| `import_parse` | `import_data_from_parse()` on the packaged Parse example |
| `assign_low_maximum`, `assign_low_thresholding`, `assign_low_mixture` | `assign_grnas()` with each method on the low-MOI example data |
| `assign_high_mixture`, `assign_high_thresholding` | `assign_grnas()` on the high-MOI example data |
| `workflow_*` (7 cases) | The full workflow on synthetic data: import, parameters, assignment, QC, calibration, power and discovery. Low and high MOI; NT and complement controls; permutations and CRT; union, Bonferroni and singleton integration; each side; skew-normal or no approximation; output amounts 1 to 3 |
| `example_low_defaults` | Low-MOI example data with package defaults (maximum assignment, NT controls, permutations), and `write_outputs_to_directory()` |
| `example_high_defaults` | High-MOI example data with package defaults (mixture assignment, complement controls, CRT) |
| `pairs_constructors` | `construct_positive_control_pairs()`, `construct_cis_pairs()` and `construct_trans_pairs()` |
| `pairs_trans_caller_mutation` | Whether `construct_trans_pairs()` changes the caller's table (known bug) |
| `pairs_cis_unsorted_positions` | `construct_cis_pairs()` with an unsorted position table (known bug) |
| `cache_reuse` | Reuse and invalidation of cached response fits |
| `cache_control_switch` | Switching from NT to complement controls after calibration (known bug) |
| `cache_control_switch_cleared` | The same with the cached fits cleared: what the fix should produce |
| `edge_duplicate_pairs`, `edge_no_passing_pairs` | Rejected requests |
| `edge_intercept_only_nt` | `~ 1` with NT controls (known bug) |
| `saved_object_0_99_0` | Loading a saved 0.99.0 object, and running discovery on it |

A further test checks that two parallel workers give the serial results
(skipped on Windows).

**Inputs.** The synthetic data keep the guide design, covariates, pairs,
removed cells and seven workflow settings of the characterization suite of the
refactor (#234). The responses are negative binomial (size 2) instead of
Poisson. With Poisson counts the dispersion estimate sits near its cap of 1000
and moves with floating-point round-off, which changed p-values between
platforms (0.052 vs 0.0051). The covariate `z` is rounded to 6 digits, so the
inputs are bit-identical on every platform. The example-data cases use the data
sets shipped in `data/` and `inst/extdata/`. The manifest records the digests
of all inputs, so a changed input is reported as such.

### Public outputs and internal state

Each observation has two parts, compared by separate tests:

- **Public** (`test-consistency-public.R`): what a user can see. This covers
  the exported getters (`get_response_matrix()`, `get_grna_matrix()`,
  `get_cell_covariates()`, `get_grna_assignments()`), `get_result()`, the
  printed summary and the files of `write_outputs_to_directory()`. It also
  covers pair constructors, exports and their arguments, errors, warnings and
  the random number generator state.
- **Internal** (`test-consistency-internal.R`): `sceptre_object` slots, stored
  under elements named `internal`. These are cached response fits, resampling
  budgets, the QC'd cell and index lists, the requested pairs with their
  pairwise QC, and the slot list.

A change that should keep results must pass the public tests. A restructuring
that moves or redefines slots, such as the move to `SingleCellExperiment`,
may run `consistency.R run --waive-internal`. That skips the internal tests,
and the pull request must say so. Adapting `consistency_observe()` to the new
object and regenerating from the reference commit restores the internal
checks.

## Known-bug cases

Some cases record behavior of the reference commit that a planned bug fix will
change. They record the reference behavior anyway, and each such case isolates
its bug, so that a fix shows up as a change to exactly the flagged fixtures.
They are listed in `consistency_known_bugs` in `harness.R` and under
`known_bug_cases` in the manifest. A failure in one of them says so.

| Case | Fix | Recorded behavior |
|---|---|---|
| `import_dense_symmetric` | A1, #4 | The lower-triangle value of a symmetric dense input is lost on import. |
| `import_cellranger_one_directory_covariates` | A2, a7b08441 | A single-directory Cell Ranger import drops the user's `batch` covariate (`lane` is kept). |
| `pairs_cis_unsorted_positions` | A2, #9 | With unsorted positions 90, 1000, 110 (midpoint 100, distance 20), only the first response is paired. |
| `pairs_trans_caller_mutation` | A2, #8 | `construct_trans_pairs(pairs_to_exclude = "pc_pairs")` turns the caller's data frame into a `data.table`. |
| `edge_intercept_only_nt` | A2, #11 | Fails with "non-conformable arguments". |
| `cache_control_switch` | A2, #12 | The discovery after switching controls reuses 4 NT-cell fits, so it differs from `cache_control_switch_cleared`. |

To land a fix that changes one of these cases:

1. Run the suite on the fix. Confirm that the only failures are the cases the
   fix is meant to change. A failure anywhere else is a regression to explain
   first.
2. Remove those cases from `consistency_known_bugs` in `harness.R`, and commit.
3. Regenerate from the fixed commit:
   `Rscript tools/consistency/consistency.R regenerate --ref HEAD`.
4. Check that the fixture diff touches only those cases' files
   (`fixtures/cases/<case>.rds`, `fixtures/tables/<case>/`) and the manifest.
   For #12, the public outputs of `cache_control_switch` should then equal
   those of `cache_control_switch_cleared`.

The other planned fixes have no case here:

- Old behavior that builds invalid matrices or reads out of bounds, which has no
  dependable output to record: #1–#3 (A1), #10 (A2), #13 and #14 (A3).
- Input validation and error reporting: #16, #17 and #27.
- Plots, and changes to the caller's own objects other than #8: #5, #7 and
  #22–#26 (A4).
- #237 (A3): no case produces a p-value above 1.

## Comparison rules

- **Doubles** are compared elementwise and relatively:
  |actual − expected| ≤ tolerance · |expected|. This holds at every magnitude,
  so a p-value of 1e-25 that becomes 1e-14 is a failure. An expected value of
  exactly 0 allows |actual| ≤ 1e-12. The positions of `NA`, `NaN` and `Inf`
  must match.
- **Everything else must be identical:** integers, logicals, strings, factors
  with their levels, names, dimensions, classes, row names, row order, error
  messages, warnings, and the R random number generator state after each case.
  Tables are stored in full (all columns, classes and row names). They are
  flattened only for the TSV copy.
- A failure lists the differing paths, for example
  `observation$value$tables$discovery_result$columns$p_value`.

### Cross-platform policy

| Tier | Relative tolerance | Used |
|---|---|---|
| `strict` | 1e-8 | On the recording platform: same `R.version.string` and platform string as `environment` in the manifest |
| `portable` | 1e-6 | Everywhere else |

`run` chooses the tier from the platform; `--tier` overrides the choice. Under
both tiers, every decision (QC, significance, pair identity, assignments) must
be identical.

- **CI.** macOS runs R 4.6.0 (the manifest's version) at `--tier strict`, and
  must pass. Ubuntu runs the portable tier as an informational job
  (`continue-on-error`) until it has passed in CI.
- **Measured on Linux.** In an Ubuntu 24.04 arm64 container with R 4.5.1, the
  suite passes at the portable tier with OpenBLAS and with reference BLAS.
- **The one strict gap.** At the strict tier, one internal value differs on
  Linux: a cached dispersion of the 0.99.0 saved object (`r5`, relative
  difference about 1.5e-8). That object holds #234's Poisson counts and cannot be
  re-recorded, because it is a file-format input. Its public outputs pass.

## Regenerating the fixtures

```sh
Rscript tools/consistency/consistency.R regenerate --ref main
```

This installs the reference commit into a temporary library and records every
case in a fresh R process, using the harness of the working tree. It then
replaces `manifest.json`, `cases/`, `tables/` and `inputs/` under `fixtures/`.
`--out DIR` writes to another directory instead.

Regenerate in three situations only:

- **The harness changed.** The manifest records the SHA256 of `harness.R`,
  `fixture-io.R` and `export.R`, and the fixture tests fail until the fixtures
  are regenerated. Regenerate from the reference commit. If the change only
  reorganizes the harness, the diff should touch only the manifest.
- **A change is meant to alter results**, such as a bug fix
  ([procedure](#known-bug-cases)). Explain each changed fixture in the pull
  request.
- **The reference moves on**, for example to a new release. Regenerate from it.

Never regenerate to make an unexplained failure go away. Regenerate on the
recording platform: with the same R version, platform and packages,
regenerating twice from one commit gives byte-identical files except for
`generated_at` in the manifest. A different R version or platform changes the
serialized `.rds` bytes, and can change the last digits of doubles.

## Fixture format

This section describes what a test in another language, such as pysceptre's,
needs to read the fixtures.

**Manifest.** `fixtures/manifest.json` contains:

- `reference`: the commit and package version that produced the fixtures.
- `environment`: the R version, platform, random number generator and package
  versions.
- `comparison`: the rules and tiers above.
- `seeds`, `inputs` (digests) and `known_bug_cases`.
- `cases`: one entry per case, with its `group`, `description`, `known_bug`
  (or null), `error` (the error message, or null), `warnings`, and `tables`.
  Each table entry gives its `file`, its `columns` with their types, `n_rows`,
  and any `omitted_columns`.
- `files`: the SHA256 of every fixture file.

**Tables.** Files are UTF-8, tab-separated, with one header row and LF line
endings.

- Each column has one of four types, given in the manifest: `integer`,
  `double`, `logical` (`TRUE`/`FALSE`) or `character`. Factors are written as
  their labels.
- Missing values are `NA`.
- Doubles are written with 17 significant digits (`%.17g`), so a correctly
  rounded parser, such as Python's `float`, recovers each value exactly.
  `Inf`, `-Inf` and `NaN` are written as such.
- Row names that are not automatic become a first column, `row_name`.
- Table names follow the observation path, joined with `__`, for example
  `warm__discovery_result.tsv`. Tables from internal state start with
  `internal__`, for example `internal__discovery_pairs_with_info.tsv`.

**Indices.** Cell and feature indices are 1-based, in the order of the input
matrices.

- `grna_assignments` lists one row per (`grna_id`, `cell`) assignment.
- `cell_covariates` has one row per cell, in input order.
- Imported matrices are `row`, `col`, `value` triplets of the nonzero entries,
  ordered by column, then row.

**Result tables.** `calibration_result`, `power_result` and `discovery_result`
are the tables returned by `get_result()`. `internal__discovery_pairs_with_info`
and `internal__positive_control_pairs_with_info` are the requested pairs with
their pairwise QC. With `output_amount = 3`, the resampled null statistics
(`z_null_*` columns) are kept in the `.rds` only. The manifest records how many
columns were omitted.

**What exists only in the `.rds`.**

- Non-tabular internal state: cached response fits, resampling budgets, the
  QC'd index lists.
- Printed summaries.
- The count matrices of the file imports. Read the original files under
  `inst/extdata/` instead.
- Trans pairs without exclusions or with targets excluded, which are full cross
  products.
- The per-cell tables of nested observations.

**Inputs.** `fixtures/inputs/synthetic/` holds the synthetic data as dense
TSV matrices (first column: the row identifier), its covariates, pairs and
removed cells. `workflow_settings.json` gives the seven workflow settings and
the arguments common to all of them. `fixtures/inputs/example/` holds the
package example data as gzipped Matrix Market files, with row and column
identifiers in separate TSV files.

**Random numbers.** Calibration pair selection and the resampling draws depend
on R's Mersenne-Twister generator, seeded with `seeds$cases` before each case.
Matching those outputs exactly requires reproducing R's random number streams.
