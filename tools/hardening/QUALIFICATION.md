# Hardening qualification, 2026-10-19

This is the preserved qualification of source **169c16e**. The subsequent
submission-diagnostic cleanup has a separate [qualification report](submission-cleanup/QUALIFICATION.md)
and archive identity; the receipts below are historical and have not been
relabelled as checks of that newer source.

All ten recommendations are dispositioned in [DISPOSITION.md](DISPOSITION.md).
The implementation is frozen at **169c16e4ba00c687c1f03fe37e9d79a11d01cd8f**.
The later evidence commit changes only excluded tooling, CI and architecture
documentation; it does not change the packaged implementation, tests or help.
No push, merge, release or hosted workflow dispatch was performed here.

## Exact final artifact

- Archive: `/private/tmp/sceptre-hardening-qualified-check/sceptre_0.99.0.tar.gz`
- SHA-256: `fb6c3ad20da63c2d2595e08b54a8916e8eed9002e50fa62e877867718afa11b4`
- Built by default `R CMD build source` from a clean `git archive` of 169c16e;
  the build successfully created the vignette.
- [Source manifest](validation/artifact-source-manifest.json) verifies all 176
  tracked R/native/help/test files byte for byte against that archive and records
  each SHA-256. [Verifier](validation/verify-artifact.py) is reproducible from the
  repository root. Final receipt edits leave all packaged source paths unchanged.
- [Build log](validation/qualified-build.log), [check log](validation/qualified-check.log),
  [00check.log](validation/qualified-00check.log),
  [installed suite](validation/qualified-installed-testthat.Rout), and
  [BiocCheck log](validation/qualified-bioccheck.log) are terminal receipts.

The older f9f8142 provisional archive and the 8cbecfb archive under
`/private/tmp/sceptre-hardening-final-check/` are historical candidates. The latter
had an unqualified `slotNames()` NOTE and advertised the rejected 3.1.7 floor.
Neither is relabeled as this final artifact; their original paths remain intact.

## Package and scientific checks

| Check | Terminal result and scope |
|---|---|
| Exact final archive, `R CMD check --no-manual --no-build-vignettes` | **Status: OK**, macOS arm64 / R4.6.0. Examples, installed tests and running vignette code pass; vignette rebuilding in check is skipped because it was built separately. PDF manual compilation is outside this invocation. Repository-index access warnings appear in the log; dependency checks still finish OK and produce no final check warning/NOTE. |
| Exact final installed suite, during R CMD check | **5,641 passes, zero failures/warnings/skips**, testthat3.3.2 and optional ondisc1.2.0. |
| Source suite at 8cbecfb checkpoint | **253 cases, 5,641 passes, zero failures/errors/warnings/skips**. [Summary](validation/full-tests-summary.json), [log](validation/source-full-tests.log) and RDS remain distinct from the exact final installed receipt. The later production change qualifies `methods::slotNames()`; the installed suite above covers it. |
| Historical matched-budget reference parity | **Seven workflows, 231 tolerance-zero comparisons**, exact results, preparation/fits, warnings and post-run RNG state. Original [source receipt](validation/parity-summary.json) is retained. The final installed-archive run has its own [summary](validation/installed-parity/parity-summary.json), [environment](validation/installed-parity/parity-environment.json) and [log](validation/installed-parity/parity.log). This compares the saved historical reference with matched budgets, not a new scientific fit of that reference. |
| Distributed package seam | **17 cases / 105 assertions**, including fresh in-memory and reopened backed R workers, direct/collect numerical and RNG equality, complete-family correction, failed/corrupt/conflicting units, ownership and runtime/build certificates. This is the R API boundary, not scheduler qualification. |
| Formatting | Air0.9.0 check passes. The mechanical 74a421a → f9f8142 proof checks all 46 R files after normalization of single-expression braces only; it does not claim raw parsed trees are identical. [Proof](validation/format-ast-proof.json) and [script](validation/verify-format.R). |
| Native source | No `src/` change from the frozen cef47da5 baseline. Native statistical algorithms and draw ordering are retained. |

## Dependency floor and support

The retained R floor is **R >= 4.1**; the qualified test dependency is
**testthat >= 3.1.8**. Both exact-source full suites use genuine, unmodified
archived testthat3.1.8 in isolated libraries:

| Runtime | Cases | Passes | Failures / errors / warnings | Explicit skips |
|---|---:|---:|---:|---:|
| R4.4.1, macOS arm64 | 253 | 5,505 | 0 / 0 / 0 | 10 |
| R4.1.3, Linux arm64 Docker | 253 | 5,505 | 0 / 0 / 0 | 10 |

Each lane reports 5,515 expectations including skips. The ten skips concern
optional ondisc/backed behavior, because ondisc is absent in those environments;
their names are recorded, not silently counted as passes. Imported-print mocking,
28 inventoried APIs/argument forms, and package/spec construction with
`R_DEFAULT_PACKAGES=NULL` pass. The ordinary user libraries were not modified.

3.1.7 exports those APIs but fails the actual inherited-print mocking required by
the sink error/interrupt tests on both R4.4 and R4.1. Those rejected-candidate logs
remain labeled separately. Its unmodified source also cannot compile under R4.6
because of removed R internals; current-R checks use current testthat. No tests
were weakened to retain an inaccurate floor. See [dependency report](dependencies/README.md),
[status](dependencies/evidence/qualification-status.json), and
[provenance](dependencies/evidence/provenance.json) for exact image/archive/source
hashes, full results and reproducible Docker/CI setup.

The standard OS jobs and pinned optional-backend job remain; R-devel, R-floor and
official Bioconductor release/devel jobs are configured. **Hosted runs were not
executed.** Local arm64 evidence does not establish Windows, x86, R-devel or
Bioconductor-container success. The coherent 2024-03-01 CRAN snapshot makes the
R4.1 dependency closure reproducible rather than resolving incompatible newest
packages. CI regenerates the API inventory before probing, and the standalone
probe regenerates a missing inventory; fresh-checkout and unknown-future-API
probes are recorded in the dependency evidence.

## BiocCheck and optional backend limits

BiocCheck1.47.26 on the final archive completed with **0 errors, 1 warning and
7 notes**. The invocation disables the R-version and Bioc-help-subscription
checks (`no-check-R-ver`, `no-check-bioc-help`); it is not an official submission
acceptance. The warning is absence of a Bioconductor infrastructure dependency.
An unused dependency was not introduced to hide it. The notes concern:

1. Literal `1:3` output-level membership checks in the execution API.
2. `cat`/`print` outside show methods, including intentional summary printing.
3. `paste` in a condition signal.
4. Three `<<-` sites.
5. Function length.
6. Missing runnable examples on the three experimental execution API help pages.
7. Lines longer than 80 characters.

The API pages document the contracts and the real regressions exercise them;
the example note remains unresolved. These are disclosed diagnostics, not a clean
BiocCheck claim.

Local R4.6 tests qualify the installed **ondisc1.2.0** interface/behavior; it has
no RemoteSha, so its version is not a verified upstream source pin. Official
upstream main/release inspection still found no supported public replacement
for the required private count operations/shared-ID constructor. The lazy gate
checks interface shape and native arities, not arbitrary future scientific
behavior. [ONDISC_CONTRACT.md](ONDISC_CONTRACT.md) records the exact inspected
upstream commit, retained semantics, source-pin recommendation and external
support/distribution gap. No new upstream-release/platform qualification or
backed CRT support is claimed. Active backing files remain immutable; external
content tamper monitoring is not introduced.

## Performance and independent review

The [bounded benchmark](benchmark/README.md) binds all 47 source/script hashes to
169c16e. Actual CRT output-level-three record ownership takes median **1.50ms**
(observed batch-average range 1.40–2.00ms), with 566,288 profiled R allocation
bytes. A separate 500-row/128-arbitrary-draw synthetic graph is explicitly labeled.
The tiny 160-source-cell/155-retained-cell, five-pair fixture has four response
and two perturbation fit calls cold and zero of each warm, with exact result,
fit and warning equality. Timings are comparable and noisy, so there is **no
speedup claim**. Rprofmem measures R allocations, not peak/native memory. This
does not qualify a large screen or arbitrarily large history.

Independent review used a clean installation of **8cbecfb** and passed twelve
focused checks verifying both reported collector fixes: conflicting valid warm
fits reject atomically, and proposals no longer repeat full cell vectors.
It also checked legitimate stale replacement and proposal ownership. See the
preserved [report](review/installed-8cbecfb/RECHECK.md) and
[verification](review/installed-8cbecfb/verification.json).
The later [narrow 169c16e review](review/narrow-169c16e/RECHECK.md) found the CI
inventory omission and confirmed the explicit namespace/floor changes. The
inventory fix is included in the evidence commit and owner-tested from a clean
temporary checkout. The reviewer did not rerun its twelve checks or certify the
final full suite/container/pipeline: this report does not upgrade that limited
review into final cross-repository signoff.

The Nextflow adapter is maintained and qualified separately. Pipeline source,
container image, scheduler, run/resume identity and direct/distributed parity
must bind this exact package archive in that separate report. Package checks
above do not establish Nextflow compatibility.
