# Internal refactor verification

Final source 0a9f9d8 is qualified. See [final acceptance](FINAL_ACCEPTANCE.md)
for the explicit source/coverage map and two bounded accepted memory costs,
and [retained evidence](results/final-0a9f9d8/README.md) for restoration.
Earlier checkpoint results below remain historical.

This directory contains developer-only characterization and qualification tooling. It is excluded from
the package tarball. `REFACTOR_SPEC.md` defines the accepted scope; live stage
results and remaining gates are summarized in `STATUS.md` here.

The source baseline is sceptre commit
`21f9ea098b69c444a884a49bc254de711968e3b5` (0.99.0). PR 233 is excluded.
The integration branch is `codex/refactor`; workers use separate branches and
worktrees. RF0 froze the original production behavior before RF1–RF9 changed
internal ownership. Public exports, defaults and the S4 layout remain unchanged.
See [ARCHITECTURE.md](../../ARCHITECTURE.md) for the resulting module boundaries.

The quantitative [maintainability scorecard](MAINTAINABILITY.md) freezes source
complexity, ownership/site inventories and repeatable change-tracing questions
before structural work. Metrics and their limitations are separate from the
compatibility gates.

## Integration checkpoints

Commit each coherent change together with its relevant tests and fixtures.
During active coding, aim for a recoverable checkpoint about every 15–30
minutes when a logical unit is ready; this is a guideline, not a timer or an
industry-standard interval. Run targeted checks before each integration commit
and broader checks at stage boundaries. Commit a stable checkpoint before a
handoff or major restructuring. Keep expectations, generators and integrity
manifests together whenever separating them would break intermediate commits.
Do not commit incomplete or unrelated work just to meet the interval. Separate
measurement/CI/documentation units when they are independently reviewable.

## Reference fixtures

Fixtures are committed under `tests/testthat/fixtures/refactor`. Tests load
saved inputs and expectations; they never regenerate them. `manifest.json`
records source hashes, dependencies, platform, RNG conventions and artifact
hashes. `api.json` is the readable export/formal/method/slot inventory; R's
generated `.__T__` method-table bindings are excluded from user API comparison.

The seven small workflows exercise both MOIs, NT/complement controls,
permutation/CRT, union/singleton/Bonferroni, all three analysis roles and output
amounts 1/2/3. They use the package's production budget selection, including
calibration's subsequent B3 update. Indexing fixtures retain noncontiguous
source cells and the distinct NT and treatment-plus-NT ordering. Other fixtures
cover warm/cold caches, changed assignments/formula/cells/side, baseline saved
objects, duplicate requests, no passing pairs and the existing intercept-only
NT-path failure. Baseline limitations are observations, not desired new rules.

The backed screen repeats the same small count pattern to 1,280 cells. That
exceeds the existing importer's 1,000-stored-entry chunk limit; it does not
establish out-of-core scalability. Scientific parity compares ordered values,
QC and results, while backed import's own RNG use and dependency warnings are
recorded separately. The selected dependency version gets exact backed
warning/RNG checks. Reattachment continues a discovery-only baseline; it is not
compared with a run whose budget was changed by calibration.

Native fixtures cover 17 prepared cases, stage thresholds, accepted and failed
skew-normal approximation, all sidedness, positive/zero B3 and selected null
segments. Only ordinary zero-based arrays are saved, never external pointers.
The test-only C++ bridge is outside production `src/` and exports no package API.
Every platform compares the full banks of exactly three active samplers with
an immutable baseline source oracle compiled in that same test environment.
The literal candidate map is Fisher-Yates -> `fisher_yates_sampler`, hybrid ->
`hybrid_fisher_iwor_sampler`, and production fast CRT -> `sample_crt_indices`.
The fast CRT comparison uses all 5,498 draws, never a diagnostic prefix.
Actual banks must come from the installed candidate namespace; expected banks
come only from the independent frozen original source. The fourth frozen
configuration, slow CRT, has no production substitute: its source/config/archive
integrity is checked on every platform and its original-oracle full bank is
compared with the historical bank on the recorded platform/toolchain. This is
archive evidence, not a fourth active candidate sampler. The controlled pointer
ownership probe covers three retained samplers x three failure phases (nine
cases). The original bridge/include/configs and historical reports stay intact. Its pinned Git-blob identity is
checked during tests. The saved platform-specific stream is an additional
check, so differing hosted compiler versions do not disable exact RNG parity.
Floating comparisons use the predeclared tolerance `1e-8`; memberships, order,
API, QC/significance decisions, stages and budgets use exact comparisons.

## Deliberate regeneration

Regeneration requires Git, a compiler, pkgload, Rcpp/BH, jsonlite, digest and the
chosen runtime dependencies. Use a fresh archive of the baseline, not a
candidate checkout. Both scripts verify baseline source blobs, reject added
source and force a clean native rebuild; they recheck source hashes after
compilation. Review the resulting fixture diff and hashes before committing.

```sh
rf0_reference_dir=$(mktemp -d)
git archive 21f9ea098b69c444a884a49bc254de711968e3b5 | tar -x -C "$rf0_reference_dir"
Rscript --vanilla tools/refactor/generate-native.R "$rf0_reference_dir" tests/testthat/fixtures/refactor
Rscript --vanilla tools/refactor/generate.R "$rf0_reference_dir"
```

Run the workflow generator last: it finalizes SHA256/MD5 digests for the complete
fixture set, including native artifacts. Package tests verify the saved manifest
and all artifact MD5 digests without additional runtime dependencies.
Git attributes preserve fixture bytes across checkout line-ending settings.

The source ondisc lane uses 1.3.5 from upstream
`589803edea2a66ad15e7f07186bee16c7739ab8c` (same source tree as the local
`4d9d75e18f21d11b6da6a54cb03d4424fcf16fe4`). Installed ondisc 1.2.0 was
separately observed to emit bounds warnings. CRAN's 1.3.5 source retains an
unbounded final-pair loop that the pinned upstream source guards; its Windows
binary emitted bounds warnings in the earlier hosted run. Version
1.3.5 alone therefore does not identify the qualified implementation. All four
standard hosted jobs install and verify the exact source pin; the no-ondisc lane
remains separate. Dependency sources must not be silently interchanged in
reference/candidate comparisons or their warnings suppressed by SCEPTRE.

## Package tests and CI

```sh
R CMD INSTALL --preclean --library=/path/to/fresh-candidate-library .
R_LIBS=/path/to/fresh-candidate-library SCEPTRE_REQUIRE_ONDISC=true \
  Rscript --vanilla -e 'library(sceptre); testthat::test_dir("tests/testthat", package="sceptre", load_package="installed", stop_on_failure=TRUE)'
```

`SCEPTRE_REQUIRE_ONDISC=true` makes a missing backed dependency a failure. The
standard check lanes pin the ondisc source revision. A separate fresh Linux lane
installs hard dependencies plus test tooling, asserts ondisc is absent and runs
the package suite. Existing Linux/macOS/Windows checks remain. The dependency
action's [documented inputs](https://github.com/r-lib/actions/blob/v2/setup-r-dependencies/action.yaml)
support the explicit hard-dependency selection. Hosted results must be checked
after publication; editing workflow YAML is not evidence of successful CI.

The hosted lanes build the verified original in a separate process with the
same numerical environment as the candidate. This accommodates pre-existing
cross-platform fit differences without changing fixtures or tolerances.
[RUNTIME_ORACLE.md](RUNTIME_ORACLE.md) documents source-byte verification,
toolchain/RNG fingerprints and debug versus release pairing. A missing or stale
reference is a CI failure. Public API, ordering, decisions and exact native
sampler checks retain their independent fixed contracts.

## Downstream stage smoke tests

`consumer-inventory.json` pins 11 pipeline scripts, their private calls, slots,
CLI arguments, artifacts and hashes, plus modules, book contracts and observed
runtime availability. It is an inspected source snapshot; execution evidence is
in separate smoke reports.

Install the reference and candidate packages into different temporary libraries.
Use identical ondisc/Arrow/R dependencies for both. Run:

```sh
Rscript --vanilla tools/refactor/pipeline-smoke.R \
  /path/to/package-library /path/to/pinned/sceptre-pipeline \
  tests/testthat/fixtures/refactor /path/to/fresh-output-directory
```

Each real script runs in a fresh R process with an isolated working directory,
explicit package resolution and verified backing-file reattachment. The report
keeps exit statuses, warnings, hashes, arguments and artifact checks; the RDS
comparison payload preserves result types/order and relevant pruned/cache state.
Compare the two `comparison.rds` files with `all.equal(..., tolerance=1e-8)` and
require exact identities/decisions separately. Standard gather correction and
nuclear raw Parquet shards are distinct contracts. The expanded harness runs
30 stages, including both nuclear trial/nontrial branches and zero-passing
shards. [CONSUMER_COVERAGE.md](CONSUMER_COVERAGE.md) records exact coverage.

`pipeline-nextflow.R` uses the same arguments to execute the complete standard
and nuclear workflows through the actual local Nextflow executor. Its explicit
qualification profile disables containers and arrays and fixes the original
unseeded pod sampler's startup seed. Historical 0.10.2 containers, the archived
default array configuration and nf-test/hosted executor profiles remain
separate qualifications. `book-smoke.R LIB BOOK_CHECKOUT OUT` runs selected
unchanged book cells without a knitr cache and preserves failed prerequisites.
Known pipeline behavior (unused nuclear threshold arguments, Parquet bytes under
the optional `.rds` suffix, discovery gather saving NULL) is preserved and
recorded rather than repaired by the refactor.

## Resource harness

```sh
python3 tools/refactor/benchmark.py /path/to/package-library \
  tests/testthat/fixtures/refactor /path/to/output --repetitions 3
```

This original compact harness fixes common BLAS/thread controls, runs one warm-up and at least
three measured fresh processes for each cold/warm case, and records elapsed
analysis time plus serial-process peak RSS. Choose dependency libraries through
the process environment, consistently for both builds. The synthetic screen is
diagnostic only. [RESOURCE_BENCHMARKS.md](RESOURCE_BENCHMARKS.md) describes the
separate actual Papalexi, difficult synthetic, sharing-sensitive CRT and backed
workloads, including sampled individual-worker/process-tree RSS and matched
clean release builds. Reports distinguish bounded resource acceptance from
unmeasured full-transcriptome or true out-of-RAM scalability.
