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
  --output /private/tmp/sceptre-EXACT-mac --ondisc present --check
```

Output must not exist. The runner resolves the commit/tree, exports the committed
source (ignoring unstaged edits), builds the vignette-containing source archive,
records its SHA256, installs with packaged tests into `output/library`, and runs
the full installed suite. `--check` additionally requires a terminal R CMD check
status without errors/warnings. Notes and test skips remain in the evidence.
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

The gate verifies the exact original archive hash and accepted cleanup archive
hash. It rejects new diagnostic categories and new practice diagnostics within
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
and backend-absent tier. `check-support.yaml` remains the PR R4.1/exact-testthat3.1.8
and Bioconductor release/devel tier, with manual dispatch and weekly refresh.
The CRAN-snapshot R-floor lane intentionally omits optional SingleCellExperiment,
SummarizedExperiment and S4Vectors, as well as ondisc. It makes no historical
Bioconductor interoperability claim. Current Bioconductor lanes resolve those
suggestions from their matching ecosystem; W05 owns present/absent roundtrips.
No existing PR check was removed. All jobs use read-only contents permissions,
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
