# Independent senior review: cleanup-plan amendment

Reviewed October 16, 2026. Outcome: **APPROVED FOR THE BOUNDED IMPLEMENTATION DESCRIBED BY THE HASH-PINNED PLAN BELOW.** The initial six material amendment requests are closed. No remaining P0, P1 or material P2 plan gap was identified. No additional user decision is needed for the adopted engineering corrections as scoped; the documented scientific/backend/duplicate-policy stop conditions remain operative.

This approves a plan, not production readiness, completed fixes, resource acceptance or pristine source. No production code was edited, tests/benchmarks run, or unsafe baseline example executed during this review. The only repository write by this reviewer is this review record. Production source still matches c7d62502757f2a17741e2ef18b9e4d006bf063dc.

## Evidence and prior-request closure

The central plan, exact guard table, complete symbol crosswalk, amended execution/core proposals, original findings, complete behavior ledger and triplet protocol were read. All 53 source entries in review-coverage.json were independently SHA256-checked against both their pinned Git blobs and the current production files: no mismatch. These are 34 handwritten R modules, 14 handwritten C++ implementations, three headers and two separately classified generated binding files. Coverage and hashes establish inspection scope/provenance, not absence of every defect.

| Initial request | Amended disposition and review |
| --- | --- |
| Post-QC group authority | Closed. CL1 validates against the current nonempty NT assignment list, not metadata rows, before combinations and derived native allocation. A group lost during QC is a required public regression. Requested/default group size is never silently clamped. |
| Exact safety semantics | Closed. Separate producers retain their safe zero domains; hybrid positive N/m prerequisites, checked N+M/M+1/products/budgets, finite selected group counts and reached-stage skew-normal minimum are explicit. Zero/nonfinite candidate rate fails informatively without a new probability floor or sampling heuristic. |
| Context/factory/intercept contradictions | Closed. One lean treatment-lookup closure is retained by purpose. Matrix-preserving intercept-only correction precedes or accompanies CL3; the old deliberate vector/error route is explicitly superseded. Saved full-source design and cache authority remain. |
| Symbol and sampler decisions | Closed. Actual callers, final names/owners and deletion paths are enumerated. No test-only compatibility aliases remain. Genuine pipeline parts, restricted/trial guide selection, pod/file helpers and result adapters retain their identified contracts. |
| Behavior disposition ledger | Closed. Every B/L finding is adopted, reproduction/characterization-gated, or explicitly deferred with its reason. Affected bug domains are distinguished from unchanged ordinary fixtures and from new methods/capabilities. |
| Ownership-copy resource evidence | Closed as a plan. Dedicated matched immutable release builds, independent inputs, correctness/nonmutation checks, meaningful triplet sizes, total/phase RSS and sorting workspace, repeated paired timings and explicit exception review are specified. The 16-byte/nnz estimate is not treated as a measured peak or passing gate. No resource exception is approved by this review. |

S3's local RAII ownership and single finalizer transfer are appropriate engineering work alongside CL1. Review the implemented success/exception ownership paths, preserving the exact existing DrawBank and external-pointer/finalizer representation. This is not permission for a new pointer class or persistence scheme.

## Exact draw-bank interpretation for implementation review

The amended guard table is approved with these explicit producer/consumer meanings, which clarify existing numerical behavior:

- Observed treatment positions are R one-based at the Rcpp entry. Their length agrees with observed n_trt, and base conversion occurs once. This equality is not a requirement on CRT row cardinality.
- A producer called with B creates exactly B rows. Production requests the checked full B1+B2+B3 bank in one call, in the existing reuse scope. Selected stage bounds must be within the bank; stage offsets and the B3=0 selected-distribution fallback stay unchanged. Do not split a full CRT request into smaller batches.
- Fisher-Yates rows each contain M distinct zero-based positions in [0,n_tot). Fixed-size inference consumes only the first n_trt entries of each selected row, requiring row.size() >= n_trt, not equality. M=0/B=0 and n_tot=M=0 retain their safe producer behavior.
- Hybrid rows each contain M entries. For every integer k with m <= k <= M, the first k entries are distinct and lie in [0,N+k). The NT-discovery caller uses that k with its aligned N+k-cell fit view. The unused suffix may contain positions outside a smaller fit view; do not reject it by scanning the entire M-entry row against that smaller view. Preserve the exact prefix/order law, N/m positive preconditions and checked size arithmetic.
- CRT rows have their actual possibly-zero size and retain their existing ordered zero-based cell positions. In use_all_cells mode, consume the complete current row; never pad, truncate to observed n_trt, remove empty draws or enforce equal row sizes. Producer probabilities/ranges are checked once, and selected-row prefix-length checks stay in the existing loop. No per-pair full-bank content scan is introduced.
- All 17 frozen native inference cases remain active unchanged, including projected multi-cell and variable CRT cardinality cases. The three active sampler mappings retain full candidate-versus-independent-original banks and the recorded-platform stored banks. The retired slow CRT configuration is explicitly oracle/archive integrity only. A missing candidate symbol is a failure; neither existence-based filtering nor routing candidate actuals through the oracle is permitted.

These contracts require tests and source readback during implementation; this planning review does not claim the guards already exist or that a new bank implementation has passed.

## Scope, migration and retained complexity

The design removes underlying duplication and uncalled modes rather than relocating them: one coordinate/view authority, direct canonical sparse counts, nonmutating owned kernels, required prepared worker inputs, explicit progress/GC effects, and deletion of unused aliases/diagnostics. Separate permutation and CRT traversal, group-specific NT fits, compact cache semantics, current mixture/tail fallbacks, saved NT representation and real pipeline adapters have concrete behavioral/consumer reasons to remain.

Adopted import ordering/modality mapping, cis distance, matrix dimensions, documented option forwarding, wrong-range variable and ownership repairs are named bounded corrections. Their new regressions must use independent intended expectations; ordinary scientific/API/RNG/schema/order fixtures stay frozen. Reproduction-gated plot repairs do not select a new significance or confidence-band method.

B5 exception narrowing and B6 removal of existing in-memory full-input densification are conditional engineering work in this pass. The actual executed exception/return-state taxonomy and aggregation order/duplicate/numerical/resource evidence are implementation prerequisites, not assumptions. Preserve numerical fallback state/arithmetic and the current aggregation output class. A conflicting numerical, duplicate, rounding or backend contract stops that dependent correction with the concrete choice; it does not authorize a new method or silent skip. New backed aggregation, Bonferroni plotting interpretation, broad malformed-input redesign and global downsampling semantics remain explicitly deferred.

CL1 precedes broad restructuring; count kernel/adapter/bindings and context/worker migrations land coherently without temporary aliases. Public API/56 slots/saved objects, supported cache/RNG/order/science and real consumer boundaries remain acceptance gates. Final immutable-source package/native/consumer/CI/resource qualification and a fresh source/deletion review remain required. Historical qualification is never relabeled as evidence for newly changed paths. This approval does not waive those gates or establish zero technical debt.

## Reviewed artifact SHA256 pins

Approval refers to the exact bytes below. Later semantic amendments need focused review; updates that merely record this outcome must preserve these reviewed target decisions.

- `REFACTOR_SPEC.md`: `5181de9c7b50380fd922de75a547f3f4ee43db603e58f2db28680aaec5de125c`
- `tools/refactor/CLEANUP_PLAN.md`: `51d5af061342e01895a656507c95fe489e9bb3e9d269ae1902b6976af1d854f6`
- `tools/refactor/cleanup/review-coverage.json`: `294d286459c2809b2782296a37a0190d14483d4218c8fd19c6e83f2b0e2f2a76`
- `tools/refactor/cleanup/core-review.md`: `79863b385f549924ee896dab02d5d74153b98cdcf15f8c76260c6b2f184502c9`
- `tools/refactor/cleanup/core-inventory.json`: `286a6604f4368d5030036c3af759d3dcc3e2e0fd8e3f7e54a3d37865f7bbef5f`
- `tools/refactor/cleanup/execution-proposal.md`: `95942e0fe2a4eec741cd2cac2c63ade0edda800228ce0a20f62762a67ed9ba4e`
- `tools/refactor/cleanup/execution-inventory.json`: `eab7742dd09a803983cdd1a672e6fe0898c6d598f095d8cc0561582cd7695854`
- `tools/refactor/cleanup/symbol-dispositions.md`: `7d382d350862ccb47112f77be091e4580421cbb16b1e5715d529d417f05a4758`
- `tools/refactor/cleanup/legacy-findings.json`: `c8a64b5d8d58ea0cfd547fac085cd158485ac661fb5f6102c45db22db8595a8e`
- `tools/refactor/cleanup/legacy-review.md`: `6773caf877ee2f200d5dc7ee3d96fcbdce62a35f08afda27ea0c73f871e2c366`
- `tools/refactor/cleanup/legacy-inventory.json`: `13173b38fdc21f03d24e1d3495d1f38943dfd8977b42e5bc5961822a12bf3d1b`
- `tools/refactor/cleanup/native-review.md`: `41601307675492521bc50c33e1e488be86b940fd01c3fa416645d9ae8aa39d95`
- `tools/refactor/cleanup/behavior-dispositions.md`: `35bc47f7ec85656102aa634114e57aacd3ec8ee362dc507d58be7ccf4b811ed0`
- `tools/refactor/cleanup/triplet-resource-plan.md`: `0a221b08edb304449712f7c5fddc587cd9de1acf1f93f7c3cb8f763b8e284697`
- `tools/refactor/cleanup/senior-review-initial.md`: `6262bc716f01070fdfdaab4db267b5c56aa901e048668b77e8476e3c1b390e5e`
