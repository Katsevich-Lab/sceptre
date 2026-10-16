# Symmetric dense conversion: reproduced baseline bug and disposition

Reviewer: rf0_senior_review. Date: 2026-10-16.

**Material P2 retained accepted-input correctness bug, separate from da6a079 copy optimization.** It exists in immutable original21f9ea098b69c444a884a49bc254de711968e3b5 (R/preprocessing_functions.R) and current da6a079 (R/storage.R). Root adopted a distinct bounded engineering correction; no new scientific-method choice is needed.

## Safe independent reproduction

The tiny runner below loads each immutable installed package in a separate R process. Independent expected response is matrix(c(1,2,2,3),2,2) with equal row/cell labels c1,c2. Labels must also permit Matrix's symmetry check; differing row/column labels can correctly avoid this path. A second fixture changes the lower off-diagonal entry to2+4*.Machine$double.eps, retaining each supplied double exactly in the expectation. Matrix1.7.5 classifies both as dsTMatrix. Two direct output orientations and the public in-memory import are exercised; all returned sparse objects are valid, indices are in range and inputs remain unchanged. No malformed old CSC or unsafe arithmetic is executed.

Both original and current lose the lower off-diagonal: expected rows(1,2)/(2,3), actual(1,2)/(0,3). Public-import response totals are1,5 instead of3,5; near-symmetric values likewise disappear. These are scientific input/count changes, not merely attributes or display. The bound libraries are original `/private/tmp/sceptre-triplet-release-hi3lt_yw/original/library` (prior matched build provenance21f9ea0) and candidate `/private/tmp/sceptre-cl2-owned-columns-pt5gsi1r/library` (build.json pins da6a079); package/DLL metadata are retained in RDS.

## Cause and required correction

Public matrix→Tsparse coercion calls .m2sparse.checking and detects symmetry with tolerance, storing one triangle. set_matrix_accessibility copies that triangle into a general dgR/dgC output without expanding it. Force general shape BEFORE sparse construction through supported public `as(as(matrix,"generalMatrix"),"TsparseMatrix")`; do not first auto-sparsify and then expand symmetry, since tolerance may already have discarded genuinely different near-symmetric lower entries. Preserve full supplied row/cell order, IDs, values, caller serialization and numeric output. Add both symmetric/near-symmetric direct and public-import regressions, while retaining rectangular/triangular/empty/integer cases and sparse/backed paths.

The simple public general conversion can add a full dense payload copy: installed .m2dense calls R_matrix_as_dense with new=1 for attributed referenced matrices, unlike old matrix_as_sparse's new=0. This is dense-input-only cost; sparse/backed input need never densify. Root chose this transparent public route plus measurement over bespoke wrapper/private Matrix maintenance. Direct as(matrix,"dgTMatrix") is explicitly deprecated on installed1.7.5. Public new(dgeMatrix,x=matrix) fails standard S4 validity because matrix does not extend the required numeric slot; bypassing validity/stripping caller attributes in place is not an acceptable copy shortcut. A custom native converter or private Matrix call would add avoidable ownership/version maintenance.

Measure the actual additional dense-copy peak/allocation/timing on matched immutable da6a079/fixed release installs; report short timers honestly and do not claim a performance pass from subsecond noise. Existing triplet ownership timings remain separate, with no changed limits or relabeled initial memory crossings. The optimized sparse ownership pin stays immutable, and the dense correctness fix lands as its own reviewed unit.

## Reproduction artifacts and source hashes

- `/private/tmp/sceptre-symmetric-dense-review/run.R`: `dfb27982436fc5798df23dd3c79fc1d8ccbf85a831b5fa02796f32de4f6cf574`
- `/private/tmp/sceptre-symmetric-dense-review/original.rds`: `4174a0018af45d23e9024028f1f79c48d6fd46e579f4e0bda79890cf6bb5cdb8`
- `/private/tmp/sceptre-symmetric-dense-review/original.log`: `72a6b9e464a8a8736e94d130ea684feb12a9905e63dd4c7b81461c543bc775be`
- `/private/tmp/sceptre-symmetric-dense-review/candidate.rds`: `4aa1a5787dceb405ef099d5cf157834a7a4df5a8eced043aaf553e3e7d376660`
- `/private/tmp/sceptre-symmetric-dense-review/candidate.log`: `72a6b9e464a8a8736e94d130ea684feb12a9905e63dd4c7b81461c543bc775be`
- `21f9ea098b69c444a884a49bc254de711968e3b5:R/preprocessing_functions.R`: `34bbbf5e65aae5246110ba53177a9208702469d3b553396ac58c520b4ac18639`
- `da6a079dc86bc18d7dd54e4ba590820eacd921ca:R/storage.R`: `3755dea9c81629b0140b03230a09d6359e57169ea02580ce4f19c1da9b4df06e`
