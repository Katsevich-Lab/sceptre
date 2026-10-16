# Final refactor evidence: 0a9f9d8

Production source: `0a9f9d879db6a52e3daef237daa1711b9d542976`.
Original authority: `21f9ea098b69c444a884a49bc254de711968e3b5`.
See the [final acceptance report](../../FINAL_ACCEPTANCE.md),
[structured acceptance](acceptance.json) and
[unchanged production/test/tool/workflow binding](production-freeze.json).
Nothing is merged or released.

## Readable entry points

- [Final source CI: six successful jobs](ci/source-ci.json).
- [Installed suite: 4,740 passing checks](installed-summary.json).
- [Actual consumer summary](consumers/acceptance-summary.json): 30 stages,
  complete standard/nuclear local Nextflow, 14/15 selected book cells with the
  inherited empty-barcode failure retained. Original outputs were verified and
  reused from the historical execution, not rerun during the final refresh.
- [Complete pre-Air resource matrix](resources/full-2ab-comparison.json) and
  [48-process actual final targeted comparison](resources/targeted-comparison.json).
  The full 160-process matrix was measured at 2ab3f5a, not 0a9f9d8.
- [Independent Air/source-carry review](reviews/air-independent.md) and
  [source proof](reviews/air-resource-source-proof.json).
- [Native score-loop cost review](reviews/score-loop-cost-review.md).
- [Separate ownership-diagnostic source review](reviews/triplet-diagnostic-source-review.md)
  and [independent interpretation](reviews/triplet-diagnostic-interpretation.md).
- [Maintainability review](maintainability/REVIEW.md) and
  [actual BiocCheck report](bioccheck/qualification-report.md).
- [Independent final-report accuracy sign-off](reviews/final-report-accuracy.md).

These are verbatim convenient copies. Private absolute paths and then-active
status sentences inside signed raw reports preserve their original context;
current completion/source authority is the final acceptance report. They are
not required restoration paths. The archive manifests map the retained raw
records to portable restored prefixes.

## Complete raw evidence

Three content-addressed archives preserve exact bytes, modes and nanosecond
modification times. Identical blobs are reused across archives; **all three are
required** for the final restore. Python 3.9+ standard library suffices for
restoration. No R, original private path, network or compiled installation is
needed to restore and inspect records. The historical xz dictionary requires
about 256 MiB decompression memory.

| Archive | Bytes | SHA256 | Retained scope |
| --- | ---: | --- | --- |
| [historical-evidence.tar.xz](historical-evidence.tar.xz) | 53,283,280 | `7f4a2eaaacf43113e472bfb9dc8926e379f0a8eab60f36b70360781cff5e055e` | 4,820 paths: original/32b/2ab builds, original oracle/checks, full resource studies, diagnostics, reviews, historical consumer/metric/Bioc evidence, failures and invalidations |
| [consumer-additive.tar.xz](consumer-additive.tar.xz) | 1,065,520 | `bc1d0db1928f44021d9c583b9a4ecf59e6a414280b44bd921e492fe032fafe97` | 1,149 final consumer paths; exact new outputs and original-reuse provenance |
| [final-additive.tar.xz](final-additive.tar.xz) | 7,239,068 | `8f96b58244fd2b1268264165021a04433f69c7c0665c251f1391f048145c36df` | 2,139 final paths including all consumer paths via blob references, all six CI logs, installed/native/full check, source/runner/build records, Air, metrics, Bioc, targeted48 and separate ownership profiles |

The [file manifest](evidence-manifest.json) binds every ordinary file here except
itself. [Historical](verification/historical-receipt.json),
[consumer](verification/consumer-receipt.json) and
[final](verification/final-receipt.json) receipts record actual archive extraction,
restoration and separate verification. Both wrong-dependency controls reject
before creating output. Final restoration independently checks all 1,149
consumer paths and nested Bioc/maintainability/CI manifests. Files are independent
copies, not hardlinks. The complete historical payload is restored separately.

From the repository root, use a fresh output directory (example):

```sh
mkdir -p /private/tmp/sceptre-refactor-evidence/stores
tar -xJf tools/refactor/results/final-0a9f9d8/historical-evidence.tar.xz -C /private/tmp/sceptre-refactor-evidence/stores
tar -xJf tools/refactor/results/final-0a9f9d8/consumer-additive.tar.xz -C /private/tmp/sceptre-refactor-evidence/stores
tar -xJf tools/refactor/results/final-0a9f9d8/final-additive.tar.xz -C /private/tmp/sceptre-refactor-evidence/stores
python3 /private/tmp/sceptre-refactor-evidence/stores/sceptre-pre-air-evidence-content/restore.py /private/tmp/sceptre-refactor-evidence/historical
python3 /private/tmp/sceptre-refactor-evidence/stores/sceptre-final-consumers-0a9f9d8-additive-content/restore.py /private/tmp/sceptre-refactor-evidence/stores/sceptre-pre-air-evidence-content /private/tmp/sceptre-refactor-evidence/consumer
python3 /private/tmp/sceptre-refactor-evidence/stores/sceptre-final-qualification-0a9f9d8-content/restore.py /private/tmp/sceptre-refactor-evidence/stores/sceptre-pre-air-evidence-content /private/tmp/sceptre-refactor-evidence/stores/sceptre-final-consumers-0a9f9d8-additive-content /private/tmp/sceptre-refactor-evidence/final
```

Repeat each Python command with `--verify-only` to recheck the existing output.
Inspect each extracted store's README/manifest/exclusions for exact path-level
scope. The final output uses `final0a/`, `consumers0a/`, `maintainability/`,
`bioc0a/`, `air/`, `air-independent/` and independent-diagnostic prefixes.
Historical output uses `final32/`, `final2ab/`, `consumers32/`, `native-cost/`,
`score-remedy/` and named review/metric/checker prefixes.

## Exclusions and limits

All raw numerical payloads/logs/profiles, failed controls, invalidated attempts
and executed runner source are retained. Installed libraries, native binaries,
disposable caches and duplicated extracted sources/full Git archives are
explicitly excluded with fingerprints and reasons, not silently omitted. Final
compiled-source text and the actual built Rpackage tarball remain. Original
source archives, exact Git pins, source transformations and source/DLL/build
fingerprints support reinstallation; reproducing runtime checks requires the
recorded tools/dependencies and scientific inputs. A timing result is specific
to its workload/environment, not a portable speed guarantee.

Papalexi cold serial memory remains a generic-gate crossing with separate user
acceptance under D2026-10-16-04. The unsafe borrowed-input diagnostic is never a
production remedy and does not establish the RSS cause. The separate earlier
8M numeric-triplet exception remains bounded. No broader exception follows.
The final partial timing attempt's false-clear receipt is preserved with an
explicit erratum; only the subsequent clean 48-process run establishes results.
Nonblocking style/tooling residue and scientific/backend deferrals are recorded
in the final report. No unqualified consumer configuration, whole-book success,
Bioconductor submission approval or universal native/OOM guarantee is implied.
