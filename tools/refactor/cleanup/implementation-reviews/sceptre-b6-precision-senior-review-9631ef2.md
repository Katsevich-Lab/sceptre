# B6 platform accumulation repair independent source review

Exact source `9631ef2b8e98a65895837b735667f8fb47f1bebb`. Read the production helper, new tests, reference fixture, additional input provenance and superseding precision amendment. No R process/build/test/benchmark during the active resource window.

**Outcome: no P0/P1/material P2 source finding in the revised mechanism/tests.** This repairs the changed summation boundary on inspection; hosted Linux and revised resource/combined installed qualification remain required. The earlier 88ed37c source/resource review and Mac study are initial historical evidence, superseded for acceptance by the real Linux compatibility finding. No old failure/fixture/tolerance is erased.

## Diagnosis and actual repair

Hosted Ubuntu x86_64/R4.6.1 produced dense sums 9007199254740994 and 0.59999999999999998 where the former Mac constants expected 9007199254740992 and 0.60000000000000009. Only the dense-input iteration failed the old test. For its unique-coordinate fixture, the original first as.matrix yields identical dense values for dense/C/R inputs, so all original inputs then use the platform's dense accumulation. The initial sparse replacement changed that boundary for sparse inputs. It was a production compatibility defect as well as an invalid universal test expectation.

The new helper restores `Matrix::colSums(as.matrix(...))` on **complete selected rows**, partitioning only the independent cell columns. A cell's entire ordered row sequence reaches one original dense accumulator; there are no rounded partial row sums, custom native/compensated arithmetic or tolerance/OS-skip workaround. The original dense input path still calls the same direct dense sum.

The existing dplyr selector, metadata repetitions, first-seen vector order, vapply numeric assembly and transpose are retained. Sparse preparation still happens once; duplicate-coordinate normalization precedes repeated metadata row selection, and logical FALSE remains zero under Matrix count semantics. The distinct storage stored-position conversion is not involved. Inputs are not sorted/mutated by this helper.

Column blocks target at most 2^20 dense elements and use `drop=FALSE`. When a selected group exceeds the target row count, one full column is the unavoidable minimum and can exceed 8 MiB; the contract honestly states that exception instead of claiming a universal cap. Selected sparse-row workspace, dense output and transpose are additional. Index arithmetic is double/wide-enough before min with n_cells and never forms a 1:0 sequence. Zero cells return shaped empty named sums; one-cell/one-row behavior retains matrices at the summation boundary. Output names originate from source column names and assignment fills existing positions.

The new helper is a focused owner of the precision/workspace contract, not an alias or a configurable framework. No new public API, S4 schema, backend capability, RNG policy or scientific estimator is introduced.

## Independent test/reference integrity

The test-only original function is independently bound to 21f9ea0. I compared its function bytes below the two provenance comment lines with `git show 21f9ea0:R/import_functs.R`: they are exactly identical, including final newline. The entire fixture SHA256 is pinned in the test. Its parent environment is baseenv and all nonbase calls are explicitly namespaced; it does not call the candidate helper. Windows text normalization is disabled for this byte-pinned source and the small RDS inputs.

The 28 captured characterization inputs have the same existing SHA256 as the earlier original evidence and reside outside RF0 goldens. Another 200 fixed-seed triplet/large/fraction/metadata-repeat inputs compare candidate with the actual original function on the current platform, preserving the caller RNG and checking input serialization. The 2^53/fraction test now compares to runtime original rather than copying Mac constants. Independent hand-specified integer/binary-fraction/logical cases remain, so expected output is not only another implementation of the candidate.

The 1025x1025 test exceeds the 2^20 target (1023 full columns plus a two-column final block), uses the rounding-sensitive values in every column, checks nonmutation and separately covers zero/one columns. It will exercise extended accumulation on hosted platforms that supply it. Exact comparison, input nonmutation, class/order and public import tests remain; no tolerance was relaxed. Frozen RF0 inputs/expected/source are unchanged.

## Qualification boundary

The earlier 532 source-evaluated checks were explicitly a mixed source-helper/older installed-public-namespace probe, not an installed 963 source claim. Its wrong-namespace initial attempt remains preserved. Root now reports fresh installed/exact replay progress separately; this reviewer did not execute those jobs or infer a hosted pass.

Repeat the original declared B6 input/repetition/safety/gate protocol against fresh matching original/candidate release builds. The old attractive Mac timing/RSS values do not qualify this changed helper, and the separately accepted numeric triplet ownership cost cannot waive its gates. A later material source change needs its own pin/readback. Final workflow/consumer/API qualification remains separate.

## Source hashes

- `R/import_functs.R` SHA256 `4a5b86dfaae73c41c4fdbf28796e7b75e83c88b226627da0dd2b69b646687fda`
- `tests/testthat/test-refactor-vector-aggregation.R` SHA256 `834c59c2e910934d461414f49ef9b00ffce5358a50f581adfee37dca7ef44b06`
- `tests/testthat/fixtures/aggregation-original-21f9ea0.R` SHA256 `9cb66ee837ae924268adc953ef8ac753d35a654f410d7aa3194af5f7d03ef6d9`
- `tests/testthat/fixtures/aggregation-inputs-6a43fa5.rds` SHA256 `f53c93eabbb30a53d5fed01788ee8ef0ef4447b94ed308540410e21baa624e2d`
- `tools/refactor/cleanup/aggregation-precision-amendment.md` SHA256 `88a84996de0a26793b6e08e392bdd2aeac0214d0450650fbc37614ba146f7f6e`
- `.gitattributes` SHA256 `23b2697c1d1dbaa94b5ac95d4a5ed5a74c5329a3b0e5d676d4993382375bb09b`
