# Compatibility refresh after Bioconductor diagnostic cleanup

2026-10-19: the unchanged pipeline implementation at
`818a156a8a53cf3e34e58b05bd37d115eb4ced49` passes the bounded compatibility checks
below with the newly cleaned package archive. No pipeline feature or statistical
implementation was changed for this refresh.

## Isolation and exact artifacts

The tested checkout is
`/Users/ekatsevi/.codex/worktrees/bioc-compatibility-refresh/sceptre-pipeline`,
branch `codex/bioc-compatibility-refresh`, created clean at `818a156`.
The unrelated uncommitted arrays/HPC/report/benchmark/demo/CI work in the original
`72cc` checkout was left untouched and is outside this qualification.

- Package source: `01bb39ccf8373ae8965a2679ee4813a63ea2d85d`.
- Package evidence HEAD: `9a4694450b5d456beedc11f3876ecb1cd84cb05c`;
  this is distinct from the packaged source revision.
- Exact owner-built archive:
  `/private/tmp/sceptre-bioc-cleaned-01bb39c/sceptre_0.99.0.tar.gz`.
- Archive SHA256:
  `606363fce32238354fc0d126d79f515973da7c8287ca43b14ce53ce2d00df4b7`.
- Local image tag: `sceptre-pipeline:bioc-cleaned-01bb39c`.
- Immutable local image ID:
  `sha256:1510c16a1d5f0da6a266f8685c4b06fa1865e0116b159956af9625dd1cc249b2`.
- Runtime: Linux ARM64, R 4.5.1, SCEPTRE 0.99.0, ondisc 1.3.1.
  Nextflow 25.10.4 / nf-schema 2.6.1 on the same Mac/Docker host as the baseline.

The build used the exact archive and reused cached dependency layers. All 176
scoped source/help/test files were independently compared byte-for-byte with
`git show` of the new source commit. Native `src/`, tests, DESCRIPTION and
NAMESPACE are unchanged from the previous 169c source. An actual network-disabled
container without host mounts verified the embedded source/archive identities
and runtime versions. Local image IDs cannot be pulled as registry digests.
The installed SCEPTRE file identity is
`8f59580aaf43aee9218bb379acda39d7ce1004e8fe65e4a2012a0a7ce98dec9b`.
The original installed-package map and session information are included in the receipts.
No image or branch was pushed.

The earlier package source `169c16e4ba00c687c1f03fe37e9d79a11d01cd8f` and image
`sha256:9d150dbba87ba7d4944b55399859ac1ce7244ec9c8ac5dd18d4afeb5ddaa1623`
remain historical checkpoints. Their receipts in [the prior qualification](validation.md)
were not relabeled as tests of this new artifact.

## Checks executed for this refresh

Both existing same-image R suites exited 0 with `REQUIRE_ONDISC=true`: nine
in-memory fixtures and two actual backed permutation fixtures each, with no
backend skips. These cover the extracted execution/collector/storage helpers:
public stage handoff, full tables and plans, continuous RNG, fits and assignments,
reversed worker/proposal order, invalid proposal rejection, atomic failure,
closed backing handles, relocation, unchanged count hashes and the input CLI.
The scientific matrix is unchanged from the prior qualification.

The actual Docker Nextflow harness also exited 0 in both representative modes:

- Reference low-NT empirical permutations: direct full-table/plan equality,
  public export, and cached resume (seven tasks cached, validation rerun).
- Distributed high-MOI CRT: direct full-table/plan equality, actual work-unit
  fan-out/collection, public export, and cached resume (17 tasks cached,
  validation rerun).

Both harnesses passed source/nested-config/seed/output-collision/unknown-parameter
rejection. Unsupported dynamic includes reject explicitly; missing/cyclic includes
fail in Nextflow parsing. A real failed worker prevents collection and completed
family publication. The literal-local external configuration policy is unchanged.

The newly installed `prepare_execution` help example executed in the new image
without network or host mounts. It genuinely invokes all three public APIs:
preparation, work-unit execution and collection, and produces a nonempty result.

A read-only comparison reused saved results from the previous 169c image, without
rerunning its scientific computation. Each image extracted its own public tables,
effective plan pairs/budgets, assignments and RNG from the nine memory fixtures
at eight saved stages. All 72 payloads are exactly identical. Package/runtime
identities and fit-cache structural identities are deliberately outside that
cross-revision comparison; within-version direct/distributed suites check fits.

The unchanged eight-case nf-test suite, schema lint and other platform
qualification were not repeated merely to duplicate the previous run. The new
actual Nextflow harness exercises the reference/distributed/export/resume boundaries.
Their previous receipts retain their original artifact identity.

## Package diagnostic qualification is separate

The package owner reports the exact new archive passed `R CMD check` with
`Status: OK`, 5,641 installed assertions, seven historical workflows / 231 exact
comparisons including warnings and RNG, and all three installed API examples.
With the original official wrapper options (`no-check-R-ver` and
`no-check-bioc-help`), BiocCheck 1.49.30 (the original submission checker) and
1.47.26 both report 0 errors, 1 warning and 3 notes, with zero new diagnostic
instances relative to
submitted source `21f9ea098b69c444a884a49bc254de711968e3b5`.
Unsuppressed local runs retain the support-site timeout and version/subscription
diagnostics; these are not claimed as clean default runs.
This pipeline task did not rerun those package checkers. The owner's original
[qualification report](evidence/bioc-refresh-20261019/package-owner-qualification.md)
is preserved separately from these pipeline/container execution receipts.

## Durable receipts and reproduction

[Execution and identity record](evidence/bioc-refresh-20261019/execution-record.json),
[scientific summaries](evidence/bioc-refresh-20261019/scientific-reports.json),
[source equivalence](evidence/bioc-refresh-20261019/source-equivalence.json),
and [original-file hashes](evidence/bioc-refresh-20261019/files-sha256.json)
bind the original build/test/example logs and resume traces in the same directory.
The saved-stage extractor is included; its temporary input/output payloads are
local diagnostic artifacts, not bundled data.

From a clean checkout of the qualified pipeline baseline (or this evidence-only
follow-up), with the pinned Nextflow/plugin installed:

```bash
container/build.sh /path/to/sceptre 01bb39ccf8373ae8965a2679ee4813a63ea2d85d sceptre-pipeline:bioc-cleaned-01bb39c \
  --package-archive /path/to/sceptre_0.99.0.tar.gz \
  --package-sha256 606363fce32238354fc0d126d79f515973da7c8287ca43b14ce53ce2d00df4b7
SCEPTRE_CONTAINER_IMAGE=sha256:1510c16a1d5f0da6a266f8685c4b06fa1865e0116b159956af9625dd1cc249b2 \
  SCEPTRE_EXECUTION_MODE=distributed SCEPTRE_TEST_CASE=high_crt \
  bash tests/run-nextflow.sh /absolute/path/to/fresh-run
```

Use the new resulting image ID if rebuilding; byte-identical image rebuilds are
not promised. Use `reference` / `low_nt_empirical_permutation` for the other
scheduler case. The harness prepares fixtures, locks and the direct reference in
that same image. Run both full R scripts inside it with the container instructions.
Existing build contexts intentionally refuse overwrite; use a fresh checkout or
preserve the previous context.

This refresh qualifies small deterministic workloads on local Linux ARM64 Docker.
It does not qualify the separate feature work, large screens, cluster/cloud,
Apptainer, AMD64/Windows, hosted CI or cross-platform numerical equality. Backed
resampling remains permutation-only. The existing legacy CLI/cache migration and
literal-local configuration restrictions remain in force.
