# Independent B1 coordinate-import source review

Reviewer: rf0_senior_review. Date: 2026-10-16.

Reviewed `226a5b4bd438a04cbe9336ce4a4e73ad85be7171` against `1eddbe8aa223e5c5cd7d70fbbe062c60bf89574f`. Source and test files matched their immutable Git blobs. Work-allocation hashes below identify the committed version; subsequent coordinator changes to that document and other work are excluded. No production edits, test reruns or benchmarks were performed here.

**Outcome: no P0, P1 or material P2 source correctness finding in this bounded implementation.** Exact immutable installed-source and resource qualification remain open. Source-evaluated checks are targeted evidence, not final package qualification.

## Source assessment

`R/import_functs.R:861-884` checks finite integral dimensions and entry count before native use, reserves capacity for the final cell pointer, validates integer/nonmissing/in-range coordinates, radix-orders owned positions, rejects repeated coordinates and builds complete numeric CSC slots. The pointer helper receives sorted cells within native integer capacity. Differences used for duplicate detection are safe within the checked coordinate range. Empty entries and zero-size dimensions are supported. Every accepted value, including explicit zero, is retained; no coalescing constructor silently sums repeats.

`R/import_functs.R:502-515` maps global feature-file rows to modality-local rows with `match`, preserving original feature-file order for interleaved modalities. Cells and values remain aligned. Parse retains its actual cell/feature file convention and creates owned zero-based columns instead of mutating an aliased vector. The shared CSC constructor does not mutate its table/columns.

Tests independently specify complete pointer/index/value slots and IDs, verify explicit zeros and formerly missing entries, serialize caller inputs, cover empty entries and range/NA/capacity rejection, exercise interleaved Cell Ranger modalities, preserve independent repeated coordinates across different batches, and reject duplicate coordinates through both owned public import routes. They check completeness as well as Matrix validity. The inspected source-evaluation recipe disables helper loading and resolves private symbols in the candidate source environment; the retained old-helper-shadowing failure is excluded. Its final CSV/log reports 46 passed with zero failure/error/warning/skip and explicitly states the native library remains the CL1 install. These logs were inspected, not rerun.

## Scope record and remaining acceptance

Duplicate rejection applies to the sceptre-owned **in-memory Cell Ranger and Parse** constructors. `R/import_functs.R:593-614` still delegates backed Cell Ranger import directly to `ondisc::create_odm_from_cellranger`. This diff neither implements nor qualifies its duplicate/interleaved-modality behavior. Explicitly record this boundary in B1 disposition/completion notes rather than claiming every file-import route now rejects repeats. Root acknowledged the amendment. A full-file in-memory preflight of that external path would defeat the backed boundary and is not required by this bounded correction; a backed-policy change needs separate scope, dependency behavior and resource evidence.

Workspace stays sparse-sized: the order vector, owned ordered cell/feature/value vectors, range/duplicate-check temporaries and output cell pointers require O(nnz + n_cells) storage with several live vectors. Cell Ranger also retains its existing full file table and per-modality subset. This establishes absence of a dense assay, not a measured peak, constant or passing timing gate. Requalify ordinary and unsorted file imports on the final immutable installed build and include the changed constructor in final resource acceptance; do not relabel historical import timings. No count-coalescing method, RNG, public formal or S4 schema change was found.

## Changed-file SHA256 at reviewed commit

- `R/import_functs.R`: `75f0bcaffd4f93e4bbff4e556d433f4fa8d5d260b81f00728283ce828b158ee2`
- `tests/testthat/test-refactor-import-coordinates.R`: `8c2aa56ea55d03327634e1c7ac94c9fbea7aec30d596396cdd416a79147c2a1b`
- `tools/refactor/cleanup/WORK_ALLOCATION.md`: `73c35e5ddd8c853699024cda99f2121b19d116dd47e952a0fbcca837bdcefac7`
