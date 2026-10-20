# Exact-source environment qualification

Use these tools in an isolated checkout. They never install dependencies and
never install the package into an ordinary library. Dependencies must already
exist; their paths and versions remain visible in the receipt. A package build
or reference-equality comparison is not statistical qualification.

## One archive, fresh installations

From the repository root, under the shared local-heavy lease when working in the
overnight program:

```sh
python3 tools/qualification/qualify.py --source . --revision HEAD \
  --output /private/tmp/sceptre-EXACT-mac --ondisc present --check --as-cran
```

Output must not exist. The runner resolves the commit/tree, exports the committed
source (ignoring unstaged edits), builds the vignette-containing source archive,
records its SHA256, installs with packaged tests into `output/library`, and runs
the full installed suite. `--check` additionally requires a terminal R CMD check
status without errors/warnings. `--as-cran` matches the standard hosted checks,
including vignette rebuilds; omit it only for the documented ordinary check lane
(`--no-build-vignettes`). Notes and test skips remain in the evidence.
Build/install/test/check logs and terminal receipt JSON survive a failed run.
Python3.9+ and R/Rscript are required. For another R runtime supply both `--r`
and `--rscript` from that runtime. Initial `.libPaths()` are captured as read-only
dependency fallbacks before redirecting R_LIBS_USER. All BLAS/OpenMP threads are
one, and compilation is capped at two jobs. Local process CPU admission still
belongs to the program's lease.

Use the exact resulting archive for the other platform lanes:

```sh
bash tools/qualification/container.sh IMAGE_ID linux/amd64 \
  /private/tmp/sceptre-EXACT-mac/sceptre_0.99.0.tar.gz \
  /private/tmp/sceptre-EXACT-amd64 present
bash tools/qualification/container.sh FLOOR_IMAGE_ID linux/arm64 \
  /private/tmp/sceptre-EXACT-mac/sceptre_0.99.0.tar.gz \
  /private/tmp/sceptre-EXACT-r41 absent
```

The driver never pulls images. Docker must already have the admitted platform
image with all dependencies. Inputs and tools mount read-only; a fresh output
mount holds the installation, archive digest, complete image identity, package
inventory, test totals/skips and terminal status. Containers have no network,
two CPUs, 6GiB memory, 512 PIDs, no capabilities, and no privilege escalation.
The image only needs Bash and R; existing qualified images lack Python. The
AMD64 lane on an Apple Silicon host is emulated, not native Linux execution.
Container driver executes installation+installed suite; it does not claim R CMD
check or independent scientific calibration. For a direct archive installation
on any host, use `qualify.py --archive ARCHIVE --output NEW --ondisc POLICY`.

## Original-submission diagnostic gate

The actual submission is commit `21f9ea098b69c444a884a49bc254de711968e3b5`.
The exact archived checker is BiocCheck1.49.30, original target Bioc3.24; the
official wrapper used `no-check-R-ver` and `no-check-bioc-help`. Original output
was zero errors, one warning, three notes. The warning is absence of a direct
Bioconductor dependency. Do not add a meaningless dependency to hide it.

```sh
R_LIBS=/private/tmp/sceptre-submission-checker-1.49.30/library \
  Rscript --vanilla tools/qualification/submission-check.R \
  ARCHIVE NEW_CHECK_OUTPUT 1.49.30 3.24
R_LIBS=/private/tmp/sceptre-submission-checker-1.49.30/library \
  python3 tools/qualification/submission-gate.py ARCHIVE \
  NEW_CHECK_OUTPUT/diagnostics-summary.json NEW_GATE_OUTPUT
```

Run from the repository root; use absolute archive/output paths. The checker
source is preserved under `tools/hardening/submission-baseline/checker-1.49.30/`;
the shown library is a pre-existing isolated local installation, not a portable
system path. Install the preserved checker into a new temporary library if
needed, preserving its version/source identity. BiocCheck may need network for
repository metadata; log environmental errors rather than suppressing them.

The checker requires a fresh output directory and binds the diagnostic-summary
bytes into its invocation. The gate verifies that summary digest, the exact
candidate archive digest, checker/target/options, and the original and accepted
cleanup reference archive hashes. It rejects new diagnostic categories and new practice diagnostics within
an inherited category, reuses the actual-checker function scanner and complete
individual long-function/line lineage gate, and compares all flagged display
source lines by file/text with multiplicity against the reviewed cleanup. That
cleanup's nine display sites already have a durable mapping to the submission.
Changed/moved expressions need explicit reviewed lineage. A lower total cannot
pay for a new long function, long line, or display line. Emitted full inventories
and CSV lineages must be retained, not just the final summary. This conservative
text gate is not a proof of unchanged function semantics.

Keep all-default checker diagnostics separately. The retained original-source
unsuppressed1.47.26 run has one environmental HTTP error, one warning and five
notes. A matched-wrapper pass must never be labelled an all-default clean check.

## Hosted tiers and integration

`check-standard.yaml` remains the PR OS/release/oldrel/devel, native-boundary,
and backend-absent tier. `check-support.yaml` retains the PR R4.1/exact-testthat3.1.8
tier, with manual dispatch and weekly refresh. Bioconductor release/devel CI jobs
were removed at Eugene's request; matched submission diagnostics remain a separate
manual qualification gate. Retained historical hosted results remain evidence for
their recorded source and environment.
The CRAN-snapshot R-floor lane intentionally omits optional SingleCellExperiment,
SummarizedExperiment and S4Vectors, as well as ondisc and rhdf5. It makes no historical
Bioconductor interoperability claim. Present-adapter qualification must resolve
those suggestions from a matching ecosystem; W05 owns present/absent roundtrips.
All jobs use read-only contents permissions,
disable persisted checkout credentials, have timeouts and cancel superseded
runs. These workflows use pull_request (never pull_request_target), no write
permissions, no external coverage token, and upload diagnostics on failure.
W01 owns the separate coverage workflow.

W10 should cherry-pick W04 commits, build the exact integrated archive, refresh
admitted container dependencies for any changed Imports/LinkingTo/Suggests, then
repeat the archive lanes plus original-checker gate. In particular, an RcppEigen
addition cannot be qualified by a pre-existing image without that dependency.
Only W10 may push the review branch/open its draft PR; PR creation triggers both
standard/support workflows. Dispatch `check-support.yaml` at that exact review
branch for an explicit repeat. Retain job URLs, terminal conclusions, artifact
hashes and actual R/platform/dependency versions. Configured jobs are unqualified
until those exact-revision jobs finish. Windows, hosted current R, native AMD64,
and Bioconductor ecosystems cannot be inferred from local Mac/arm64 receipts.

`check-scientific.yaml` is the separate dispatch-only heavy evidence tier. Choose
one fixed W03 scenario and40/200 independent datasets; it installs a hashed source
archive into a fresh library, uses600cells and serial inference, and uploads the
DGP manifest, raw outcomes, failures and Monte Carlo intervals. It does not run
on PRs or award a statistical-validity badge. Execution success means the harness
finished; a40-replicate pilot cannot certify nominal calibration. Deliberate
misspecification and tied/sparse stress cases retain their separate interpretation.
Its archive is built without vignette execution and is not a release gate.

Hosted baseline release/oldrel --as-cran checks exposed a genuine undeclared test
import: the backing corruption regression calls rhdf5::h5read/h5write. rhdf5 is now
a meaningful Suggests dependency; future backed test environments should install
it. This can legitimately remove the old no-Bioconductor-dependency warning.
The submission gate permits diagnostic removal using a message multiset subset,
while retaining strict bounds on all remaining flagged source instances. The
CRAN-snapshot floor intentionally omits this optional Bioconductor test backend.
The data-raw generator remains in Git but is excluded from the package archive;
packaged data files and runtime code are unchanged. Actual hosted repair requires
an exact new-revision --as-cran run; the old failing jobs remain historical facts.

Gate-policy regressions are reproducible without compilation using
`R_LIBS=CHECKER_LIBRARY python3 tools/qualification/check-submission-gate.py ARCHIVE SUMMARY NEW_OUTPUT`.
The seven fixture cases include accepted warning removal and rejected new or
duplicated messages, archive mismatch, stale summary and weakened options. Their
mutated summaries are synthetic policy fixtures and must not be used for admission.
Optional SeuratObject from W05's subsequent adapter scope is deliberately absent
in the historical CRAN floor lane; current present/absent adapter qualification is
separate. The installed-suite receipt records its availability on future runs.

Validate workflow expressions and context availability with `actionlint` before
publishing CI edits, not just a generic YAML parser. The existing isolated local
binary is `/private/tmp/sceptre-ci-actionlint-1.7.12/actionlint`; invoke it against
all `.github/workflows/*.yaml`, including coverage. Job-level env cannot use the
runner context. Resolve dynamic runner paths in a step and persist them with
GITHUB_ENV. The PR240 zero-job support failure and repaired exact-head staged
validation are retained separately from package/runtime qualifications.

The floor CI artifact uses an explicit whitelist of terminal JSON, complete
build/install/test/check logs and environment/version tables. The prior entire
output-directory upload included private installations, duplicated build sources
and opaque RDS and produced a110MB artifact. Archive identity remains in the
receipt; full local output is unchanged. Artifact minimization does not remove a
test, diagnostic or qualification step.

DelayedArray is likewise a genuine optional test dependency: the SCE rejection
fixture constructs a delayed assay to verify that backed input is refused. Its
new direct Suggests declaration belongs to the final package owner. The CRAN-only
historical R-floor lane omits it with the other optional Bioconductor backends;
current Bioconductor/present-adapter checks must retain and exercise it. This
avoids requesting a Bioconductor package from the CRAN snapshot and makes no
claim that the DelayedArray fixture ran in the absence lane.
