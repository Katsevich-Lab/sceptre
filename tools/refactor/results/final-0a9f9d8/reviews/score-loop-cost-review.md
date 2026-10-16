# Independent bounded native-loop remedy source review

Status: **final signed source review; no P0/P1/material P2 finding**. Base is immutable `32b2076b8ce6fd27cfa1a6326890342c241dbd27`. Independently reviewed the one-file patch and bound the final immutable worker commit `de27568195882d456f692e4be1b9e42f67f11dc0` and integration commit `2ab3f5a2cf3d9f19242cbe8c419ada913aa721a8` below, without editing production or running R/compiler/test workloads.

**No P0/P1/material P2 source finding identified in the inspected patch.** The change is a bounded representation-access repair, with no new method, numerical policy or extension layer. Approval is confined to the exact reviewed two-expression access change and explanatory comment; renewed gates must qualify the committed integration build.

## Exact change

Only `src/inference_scores.cpp`, within `null_score_segment()`:

1. Acquire `const double* D_data = D.begin();` before allocating the native segment result.
2. Replace the existing `D(i, curr_vect[j])` read with `D_data[i + static_cast<R_xlen_t>(D_nrow) * curr_vect[j]]`.

All explicit shape/cardinality/budget/segment/prefix guards, stage offsets, bank sharing, loops, accumulation order, floating arithmetic, `a`/`w` accesses and ownership arrangement are retained. No global bounds-check flag is disabled. No unrelated function, header, binding, test or scientific fixture change is needed by this source remedy.

## Pointer, shape and coordinate proof

- `D` is the existing const `NumericMatrix&`, alive/rooted by the surrounding Rcpp input/caller throughout evaluation. No code resizes or mutates it. The pointer is borrowed read-only; the change neither clones nor persists it. C++ segment allocation cannot relocate the existing R payload.
- The normal owning entry `run_association_test()` checks aligned `a/w/D` columns, aligned `y/mu`, observed cardinality/positions and budget arithmetic before result ownership. Its local observed-position and XPtr wrappers finish before tail creation. The staging caller keeps the same checked arrays for every segment. `D.nrow()` reads Rcpp's cached row count and is unchanged.
- For `0 <= i < D_nrow` and a consumed fit column `c`, the new offset is exactly Rcpp Matrix's own column-major formula, including the **same R_xlen_t widening before multiplication**. It is not signed-int `nrow*c` followed by a cast. A valid matrix and valid column put every accessed offset within the existing matrix payload. Zero rows skip the D-read loop; empty treatment rows skip its inner loop.
- Ordinary/complement permutation banks use Fisher-Yates positions in `[0,n_fit)`, with the fixed observed treatment prefix. NT calibration's Fisher bank uses the ordered NT-only universe. The unchanged sampler guards establish representable sizes.
- Hybrid permutation banks deliberately store a largest-size row whose suffix can belong to a larger universe. Only the size-k prefix is consumed for the current N+k treatment-plus-NT view. The patch keeps the prefix-length guard and loops unchanged; it does **not** scan or reject a legitimate suffix against a smaller view.
- CRT banks append the current fit-cell position as each cell is processed. Their complete variable-cardinality rows therefore contain valid current fit columns, including empty rows. The entire-row cardinality and integer-capacity guards remain unchanged.
- R workers obtain D from the same ordered response/design/fit view used to build these banks. The raw access does not alter source/analysis/NT/fit coordinate composition or draw-bank construction.

Caller inventory was checked directly: production `null_score_segment()` is called only by `evaluate_staged_tail()` for its three existing segments; production staged evaluation is called by the checked association entry. The developer component bridge has direct `component_segment()` and `component_tail_record()` paths, both with the existing shape guard and supplied/frozen valid-coordinate banks. Package native fixtures also use `refactor_native_pointer()`, which validates zero-based range/integer type/uniqueness before constructing its test bank. The bridge itself is not a public arbitrary-bank validator. A malformed fabricated private pointer/draw list is not made a supported input by this change.

The prior indexed access's incidental Rcpp cache check is not the provenance for supported bank bounds: producers/caller contracts are. This patch removes only that hot read path, retaining the existing explicit guards. No new validation or changed malformed-private-input policy is proposed.

## Supported minimum accessor behavior

I independently read the retained official-tag Rcpp 1.0.9 headers and current installed-header copies. The const accessor chain is:

- 1.0.9 `Matrix.h:112` -> `Vector.h:336` -> `traits.h:43`.
- Current `Matrix.h:112` -> `Vector.h:333` -> `traits.h:45`.

In each version the numeric matrix's const `begin()` delegates to its Vector const `begin()`, then returns the cached `start` pointer. It does not call R, allocate, coerce, preserve, or materialize ALTREP at access time. Cache initialization/materialization belongs to the already-existing Rcpp construction boundary. The minimum Matrix offset at `Matrix.h:174` is exactly the new formula. The patch therefore introduces no new R/materialization boundary under the surrounding tail owner or local segment owner, and no new R/Rcpp minimum requirement.

Exact inspected header SHA256:

| Header | SHA256 |
| --- | --- |
| Rcpp 1.0.9 Matrix.h | `b9f48c70a8cdfd8370b6879b90fb05d24f99302f4da8fb6fee1a393e7b130ce5` |
| Rcpp 1.0.9 Vector.h | `faaa1c4fb74bbf54efa27d5546ddc83bf6a1f690889c90a0c1ac4bf4dbd65f17` |
| Rcpp 1.0.9 traits.h | `c608f495889162e4e7eea79fab6f0040f4f34cfb12c1692fec09e4989d27ce1c` |
| Current Matrix.h | `b9f48c70a8cdfd8370b6879b90fb05d24f99302f4da8fb6fee1a393e7b130ce5` |
| Current Vector.h | `9da784c28970aefee8e0f8904486b4f2b9ce2e1bff48d08047ab7c9aef5011ca` |
| Current traits.h | `cdaacbb8f195766b44f5ecf4f8b6ed19bfd0b49e2c81444ba4ba1e254490d765` |

The files live under `/private/tmp/sceptre-native-cost-32b2076/`; `evidence-manifest.json` records the matching hashes and official `RcppCore/Rcpp` 1.0.9 URLs. I verified those hashes, not merely the filenames. This is source support verification, not execution on minimum Rcpp.

## Numerical/ownership judgment and limits

The pointer reads the identical D elements in the identical sequence, preserving every addition/multiply/sqrt and result order. No sampler, threshold, tail fallback or fit/cache stage changes. Existing F2 result continuation/raw-leaf cleanup and native-only scopes remain intact; generic pre-leaf/generated Rcpp precious preservation and universal allocation safety remain outside that earlier proof. The pointer adds no heap ownership or R cleanup obligation.

The native isolation report `/private/tmp/sceptre-native-cost-32b2076/REVIEW.md` explains the required remedy: exact installed replay reproduces the 32b2076 cost, and same-kernel access isolation finds the raw-D variant sufficient while retaining checks. I read its scope/provenance. I did not rerun the native replay or microbenchmark, and do not infer which compiler optimization caused the difference. **The microbenchmark is not whole-workflow CPU/RSS qualification, and the previous 32b2076 resource failures remain failures.**

The committed-source readback is complete. The required affected frozen/native/component/lifetime/workflow/resource/platform gates must qualify the explicit new immutable integration build. Whole-source reviews at cc01683/32b2076 and their deferred P3 follow-up remain preserved; this report neither overwrites them nor declares the refactor ready.

## Committed-source binding

Read back directly with Git on 2026-10-16:

- Worker commit: `de27568195882d456f692e4be1b9e42f67f11dc0`.
- Integration commit: `2ab3f5a2cf3d9f19242cbe8c419ada913aa721a8`.
- Both commits have identical `src/inference_scores.cpp` Git blob `dd45952e22cf22c44823e552c478d69b4dafb794`.
- Integration-vs-base changed-path inventory contains **only** `src/inference_scores.cpp`; the exact worker diff is 3 insertions/1 deletion and matches the reviewed patch above.
- Final source SHA256 is recorded in the companion binding JSON. No source/header/test/binding expansion is approved by this review.

Signed by `/root/whole_package_review`, independent source reviewer. This closes the bounded source-review request. Runtime gate results reported by the implementer (829 installed checks, 112 component checks and 959 controlled lifetime cases plus the separate negative control) were not rerun by this reviewer. Root's new integration archive/build and full resource matrix remain separate acceptance authorities. I will hold runtime/compiler/profiling/broad hashing/archive work during that matrix.
