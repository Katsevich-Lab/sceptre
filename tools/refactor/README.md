# RF0 characterization checkpoint

This directory contains developer-only baseline tooling. It is excluded from
the package tarball. `REFACTOR_SPEC.md` defines the accepted scope; live stage
results and remaining gates are summarized in `STATUS.md` here.

The source baseline is sceptre commit
`21f9ea098b69c444a884a49bc254de711968e3b5` (0.99.0). PR 233 is excluded.
The integration branch is `codex/refactor`; workers use separate branches and
worktrees. Production R/C++ code, exports, defaults and the S4 layout are
unchanged in RF0. The first structural unit follows the RF0 review checkpoint.

The quantitative [maintainability scorecard](MAINTAINABILITY.md) freezes source
complexity, ownership/site inventories and repeatable change-tracing questions
before structural work. Metrics and their limitations are separate from the
compatibility gates.

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
Every platform compares full sampler banks with an immutable baseline source
oracle compiled in that same test environment. Its pinned Git-blob identity is
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
separately observed to emit bounds warnings; dependency versions must not be
silently interchanged in reference/candidate comparisons.

## Package tests and CI

```sh
SCEPTRE_REQUIRE_ONDISC=true Rscript --vanilla -e 'devtools::test(stop_on_failure=TRUE)'
```

`SCEPTRE_REQUIRE_ONDISC=true` makes a missing backed dependency a failure. The
macOS check lane pins the ondisc source revision. A separate fresh Linux lane
installs hard dependencies plus test tooling, asserts ondisc is absent and runs
the package suite. Existing Linux/macOS/Windows checks remain. The dependency
action's [documented inputs](https://github.com/r-lib/actions/blob/v2/setup-r-dependencies/action.yaml)
support the explicit hard-dependency selection. Hosted results must be checked
after publication; editing workflow YAML is not evidence of successful CI.

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
nuclear raw Parquet shards are distinct contracts. Nuclear smoke currently uses
the trial branch with fixed pod maps; it does not cover every nuclear branch.

These are stage tests. The pipeline's current test profiles contain NULL
placeholder inputs, Docker has no running daemon locally, and nf-test is absent.
Full Nextflow and historical 0.10.2-container execution remain separate gates.
Known pipeline behavior (unused nuclear threshold arguments, Parquet bytes under
the optional `.rds` suffix, discovery gather saving NULL) is preserved and
recorded rather than repaired in RF0.

## Resource harness

```sh
python3 tools/refactor/benchmark.py /path/to/package-library \
  tests/testthat/fixtures/refactor /path/to/output --repetitions 3
```

The harness fixes common BLAS/thread controls, runs one warm-up and at least
three measured fresh processes for each cold/warm case, and records elapsed
analysis time plus serial-process peak RSS. Choose dependency libraries through
the process environment, consistently for both builds. The synthetic screen is
diagnostic only. Realistic difficult-fit data, a backed workload exceeding an
explicit dense-memory budget and fork-worker/process-tree RSS are still required
before making resource-regression or scalability claims.
