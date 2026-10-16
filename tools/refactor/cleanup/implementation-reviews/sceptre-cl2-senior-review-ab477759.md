# Independent CL2 source review

Reviewer: rf0_senior_review. Date: 2026-10-16.

**Outcome: no remaining P0, P1 or material P2 source finding in this bounded CL2 change.** This is source/design review, not resource or release acceptance. The required matched triplet ownership resource experiment remains open; no correctness-cost exception is approved here.

## Source and scope

Reviewed commit `ab4777591cb31906f1c154c1825f0a6bed8c7c6a` against parent `489ca7c` in `/Users/ekatsevi/.codex/worktrees/refactor-ci-oracle/sceptre`. HEAD was verified and every added/modified file below matched its immutable Git blob. Deleted implementations were reviewed from the parent. Frozen test fixture paths have no diff. The five untracked triplet-resource harness files are ongoing worker work and are excluded from this committed-source verdict. No production edits, test reruns or benchmarks were performed by this reviewer. The worker's 1,059 focused checks are reported evidence, not checks independently executed here.

## Findings and rationale

1. **Canonical coordinates and retained counts are coherent.** `src/sparse_counts.cpp:65-116` consumes one-based analysis group positions and an ordered analysis-to-source cell map, converting owned C++ copies once. The private CSR loader maps source presence to every analysis position, retaining repeated source cells. Stored zeros/fractions count by stored position, duplicate stored coordinates collapse to presence as before, repeated group membership contributes repeatedly, NT control totals sum group counts rather than unioning overlaps, and repeated requested pairs are retained. Response sorting is checked at the private native boundary; `R/storage.R:104-146` restores caller request order. There is no per-pair NT-local round trip or full-bank scan.
2. **Native ownership is explicit and counts remain bounded.** Rcpp converts nested groups, request indices, source cells and mitochondrial rows into owned `std::vector` values. No R slot/group vector is decremented. Coordinate bounds are checked before dereferencing these owned positions; equal pair-request lengths and sorted responses are checked. Group membership sizes bound individual integer counts, and the NT sum accumulates in int64 before checked integer narrowing. Empty groups/views/requests retain defined integer zeros and matrix dimensions. CSR slot validity is an explicit Matrix-input precondition, not a new general raw-slot validation API.
3. **The backed adapter retains the necessary compatibility boundary.** `R/storage.R:107-129,161-169` keeps disposable R group copies and analysis-to-NT-local conversion only for the unchanged private ondisc APIs that require them. The sparse branch receives direct analysis groups. Nuclear/pair-only summary placeholders and output column selection are unchanged. The actual adapter/caller migration and generated registration are consistent, and superseded native symbols/implementations are deleted rather than retained as aliases.
4. **Triplet ownership correction preserves the conversion result.** `set_matrix_accessibility` builds an owned data.table before sorting instead of passing caller-owned triplet slots to destructive `setDT`/ordering. Orientation, stable duplicate ordering, stored zeros/fractions/logical FALSE handling, pointer construction, dimnames and output sparse classes remain as before. No dense realization is introduced. Compressed orientation fast paths retain their previous behavior. `sparse_pointer_to_indices` is a coordinated private rename with the same implementation, not a new conversion algorithm.
5. **Cell summaries preserve their distinct value semantics.** `src/sparse_matrix_functs.cpp:39-95` now owns mitochondrial positions; active indices are range checked and converted once. The loaded-value >0.5 count, value sum, repeated mitochondrial rows, first maximum tie, disabled zero summaries and zero-total ratios are unchanged. The count-presence kernel is not incorrectly reused for these value summaries.
6. **Tests enforce the advertised invariants independently.** `test-refactor-storage-ownership.R` computes expected stored-coordinate counts independently and checks overlapping NT groups, repeated source/group/request positions, empty/single/repeated response selections, input serialization, mitochondrial repeat calls and exact sparse slot/alias preservation for numeric/logical/repeated-coordinate triplets in both orientations. Existing storage contracts cover real backed/reopened stores and retained ondisc copying, while generated native bindings and renamed calls are updated coherently. These are relevant contracts rather than only candidate-versus-candidate comparisons.

## Qualification boundaries and complexity

The fused count kernels retain the all-response traversal and output shape, source/analysis Boolean scratch, owned coordinate-vector copies and group-membership traversal. They remove the redundant sparse NT coordinate round trip without adding per-pair mapping/scanning. Their memory includes the existing all-response group result and sparse-coordinate/group inputs; no claim of smaller asymptotic output storage is made.

The triplet correction deliberately adds owned O(nnz) sparse column storage plus sorting workspace. Its total peak RSS, allocations and timing are not established by source inspection, a 16-byte payload estimate, the focused checks, or historical resource reports. Complete `cleanup/triplet-resource-plan.md` on exact immutable release builds and independent original/candidate inputs, then apply its gates or record an explicitly reviewed bounded correctness-cost exception. General sparse/backed/fork resource and final installed/consumer/CI gates also remain the coordinator's acceptance work.

This review does not expand the valid-Matrix precondition, promise arbitrary malformed raw-slot safety, add backed capabilities, or certify all pre-existing conversion edges. No scientific method/default, RNG draw, public API, S4 schema or frozen fixture change was found in the bounded diff.

## Reviewed plan hashes

- `tools/refactor/CLEANUP_PLAN.md`: `fae06927ca2ec4ae0c49a50a858f61b791ee323320471d049eeede9c6160c7b9`
- `tools/refactor/cleanup/symbol-dispositions.md`: `7d382d350862ccb47112f77be091e4580421cbb16b1e5715d529d417f05a4758`
- `tools/refactor/cleanup/senior-review-amendment.md`: `614bbc3e617421cb4df98ae47fd1b1ed49f3db3ab2493adefea7672d2fa36476`
- `tools/refactor/cleanup/triplet-resource-plan.md`: `0a221b08edb304449712f7c5fddc587cd9de1acf1f93f7c3cb8f763b8e284697`

## Immutable changed-source hashes

| File | Change | SHA256 at reviewed revision |
|---|---|---|
| `R/RcppExports.R` | M | `38941dba2ff7500d8dd5d4991bf30d4c0ededac8a1dbd9d6c98ccfc7a2c40e75` |
| `R/storage.R` | M | `2ae00632454136c7bec4ee4afbdc34ddd1796f4383546c712ce9854902986593` |
| `src/RcppExports.cpp` | M | `7e4484c794859a03e7624316bde3a71af7c8a518f6dad3b50c6eb67fd19e2f5b` |
| `src/negative_control_functions.cpp` | M | `7991bacbb68a6b607485e65539505bf3ca04219869c6f9f46804c553a3c63bcb` |
| `src/pairwise_qc.cpp` | deleted | parent SHA256 `bab7f126aa5380d92c0a8ef3360702eba2d1e0ee159e48b2ccc2664067da9b01` |
| `src/shared_low_level_functions.cpp` | deleted | parent SHA256 `f5bb2f5d2d782f0d60de91c7ae8c2bb78b2e81ac30c07b6695a8a583c345d5e5` |
| `src/shared_low_level_functions.h` | M | `1c93487d71a73a11698b5a1a15ad8a16848c62456ed767511462fc3dfa31f879` |
| `src/sparse_counts.cpp` | A | `3f78bd9412b1fe466c98240546814a25319d1470f61ed88ec504ff2224ef21d1` |
| `src/sparse_matrix_functs.cpp` | M | `7168d827429791dd340434435b8895f21fd70e197e401e4cbf21c27aefe52f04` |
| `tests/testthat/helper-refactor-storage.R` | A | `65e73c63d4a87bf668bcd7258f5d32ea9db61186b4d9f10a20b2b7efa4dc45da` |
| `tests/testthat/test-refactor-indexing.R` | M | `c61becc1b14f3066c15d88431b9cc5b11c4021fd320626f5525d9d84cd1dbbfa` |
| `tests/testthat/test-refactor-storage-ownership.R` | A | `cc04b357744a3b9354e943ce2905200ecbb96c6e7ccd63b5ca7002a0c3802e23` |
| `tests/testthat/test-refactor-storage.R` | M | `b796b6735793b5184e8159d47c29a5936734ab88dfa4dd98e00dd85def4a5c9c` |
