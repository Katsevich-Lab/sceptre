# Independent native output lifetime source review

Production pin:11502f871ebeb78d7cf084c45ffdb895d74ac937 versus dc093d27acefd2b3897eaf8b5a0bb4e9cbe7282f. Includes30e7aef, e65a6ed and11502f8. Read-only source/design/probe review; no R/compiler execution by this reviewer. Source checkout clean at review. This report qualifies only the specified native payload/scratch lifetimes and R conversion leaf, not universal Rcpp/OOM safety or final package acceptance.

## Source outcome

Source-approved for the bounded repair. No remaining P0/P1/material P2 found in the production diff. Two original allocation-longjmp findings are addressed in source: calibration sampling result allocation and selected association-tail result conversion. Exact frozen/installed/component/runtime qualification remains separate. The probe evidence retention finding below remains open and is not hidden by this source approval.

The repair uses a small shared continuation guard with a precise obligation: protect a public R continuation before native ownership, balance it on C++ compute failure, transfer the still-protected token without an R operation, end the complete compute scope normally, then enter R_UnwindProtect. Minimum R4.1 and current implementation install the preprovided-token cleanup/SETJMP before the raw leaf. The normal result is rooted in the continuation through cleanup; final UNPROTECT immediately precedes raw SEXP return. Rcpp's generated RObject assignment protects its incoming wrapped value. No native payload survives that outer preservation boundary.

Groups allocate/materialize the exact matrix before unordered_set/sample/copied-row owners. Unchanged draw/rejection/insertion/iteration loops write through a cached pointer and end before return. Allocating the matrix earlier changes bounded failure ordering but introduces no new random draw or numerical policy.

Pairs materialize input pointers before the native scope. Their unique_ptr owns exactly the existing grown accepted-result vectors; the permutation scratch has the same extent and loop order. The released payload enters a raw-C-API exact-length list conversion. Its no-throw/allocation-free cleanup deletes it exactly once on normal or R-jump return. No request-sized four-vector R capacity is allocated for all-rejected requests. Counts/thresholds/partial Fisher-Yates stream, accepted ordering and one-based output coordinates remain the same.

Inference computes effect/score and ends observed-position/XPtr temporary owners before owning TailResult. Numeric Rcpp input caches are materialized at parameter construction and read cached pointers. The owned tail retains only the actual selected segment, preserving fitter sorting/reflection and B3=0 behavior, not unreached budgets or the bank. The exact raw leaf preserves p/z_orig/fc/se as double scalars, integer stage, three SN values and optional ordered resampling_dist, with unchanged names/order. No local precious/native owner or C++-throwing operation exists in either raw leaf. Cleanup frees the native payload before R continues its jump or the generated wrapper preserves the successful R result. Caller parameters/borrowed Rcpp wrappers are expressly outside the stronger cleanup guarantee.

Review/probe findings corrected in this pin:

- e65a6ed widens raw column strides to R_xlen_t before multiplication, retaining prior Rcpp Matrix::offset support for valid long matrices. This avoids new signed-int overflow without allocating a huge test matrix.
-11502f8 removes duplicate null_score_segment shape validation that called D.ncol()->Rcpp dims Shield/PROTECT under a live native tail. Sole production entry validates before ownership; both direct component bridge entry points now do so too. Contract explicitly places that precondition at the owner. D.nrow() is a cached scalar in supported Rcpp; element reads use preinitialized numeric caches. This is coherent validation ownership, not a second fallback.

Native arithmetic in evaluate_staged_tail/inference_tail is unchanged. Observed score/effect expressions and permutation/rejection loops retain summation/RNG/order. Private pair return declaration changes List to SEXP only to avoid an additional result owner; generated .Call names/arity and R wrappers/schema remain unchanged. The helper is small and shared by the actual three conversion callers; no general allocation framework or eager max-budget buffer was introduced.

## Probe assessment and outstanding evidence finding

The actual-source interposition is suitable for this scoped proof: only vector/unordered_set allocators, copy/move-aware payload lifetime markers and absolute include paths are changed; arithmetic/control flow/production owner scopes stay present. It monitors R alloc/coercion/materialization/names/attrs/protection/preservation/callable boundaries, native allocations and native payload/scratch live counts. Raw errors occur only in unowned or cleanup-protected states; an unsafe-order negative control C++-fails before a deliberately leaking jump. Success workloads cover accepted/all-rejected pairs and selected-stage tails, including huge unreached B3 and compact output. The driver checks input/RNG stability and exact schema, with numeric tolerance stated1e-8. This is not an exhaustive allocation taxonomy or a scientific numerical proof.

Worker's passing log reports959 controlled cases across12 workloads; no actual OOM. The log and RDS remain preserved, as do two failed earlier attempts (probe accounting and actual duplicate dimension-protection boundary). No unchanged run is relabeled as a new source qualification.

Material P2 EVIDENCE RETENTION remains open at this pin: test-result-ownership.R writes six transformed source files under output basenames, then final file.copy(paths, output, overwrite=TRUE) overwrites them with originals. I verified all14 final archived originals match11502f8, which confirms the transformed files are missing rather than preserved. The RDS retains their pre-overwrite MD5s. Required remedy: separate original/transformed directories, restore exact transforms from immutable source and original include paths, match every recovered byte digest to recorded pre-overwrite hashes, retain recovery manifest, and fix/pin the reusable driver. A matching recovery does not require discarding/repeating valid computational evidence; any mismatch requires investigation/rerun. This is a validation artifact defect, not a discovered native numeric failure.

Root also requests a source-verified original oracle boundary for the standalone driver. That harness change and its bounded qualification must be reviewed separately before full probe evidence closure.

## Source/evidence pins

- `src/inference.h`: `25f6d04ed74af62f3a01c1240e95762306198e738f06750a6bbd4a37124abd78`
- `src/inference_scores.cpp`: `0627bfc90af01f573fa6f43922eb81ba3b2c08dc90bdf99e660fb8f65c3091f1`
- `src/inference_tail.cpp`: `7f5466d29de377966c8e0ff3b37a300aa2c05a1c345b82271d32c00cbdb65102`
- `src/inference_stages.cpp`: `1ff3732af360b79aeabf34e2ea10f3b336304b894a658618992152227d8d2b00`
- `src/negative_control_functions.cpp`: `5b21c748f9d45142c14ed1ecbc0a64953774daa6fde7014a06ada03483b56539`
- `src/low_level_full_test.cpp`: `8f9632f650f92662b774a9a5e325625d069fb675a9220e20b974c06742d9064a`
- `src/result_lifetime.h`: `f63b2ad2b6d4bf8c03f926c4c7020eb29362b7298d07434d2b97be1cc68dfe3a`
- `src/RcppExports.cpp`: `e471a595c9fc39e34e1d15c6aac7769eeabf8f6e1293541e1af703c894455adb`
- `tools/refactor/native-components.cpp`: `8e864c19b17d301bab5fcd369a99c6724b464493a6f252ae9aec050de7165d06`
- `tools/refactor/result-ownership.cpp`: `1dc834350de04dba62d03150f7ffbf7954248f5356769e89440d43e46bd5e364`
- `tools/refactor/test-result-ownership.R`: `91390e23ff88639a6e8d2fdaebd20e2c6f3b7cbbf5fa9cf63153b82618aae26e`
- `tools/refactor/test-native-components.R`: `7cdf8a36ea0e254362295f1354155bbd285e4e1c9962344bcf37e2afc5ca0dac`
- `/private/tmp/sceptre-native-lifetime-qualification-probe.log`: `20528ece628fdad3c87f8703c0f3832a4fbd3ca36652c05ccb51a533fc206b13`
- `/private/tmp/sceptre-native-lifetime-qualification-probe/result-ownership-results.rds`: `149c2fefeb7799c797d08c501fa438d28bc454c7c597dfd8a9d9132f7668971e`
- `/private/tmp/sceptre-native-lifetime-source-checks-11502f8.json`: `1e77197e136dfe572c3beecaf8aab7a9db4e51487b66be2a670b269f306f55ff`
- `/private/tmp/sceptre-native-output-design-senior-review.md`: `8b5dad334e8091e054e230445f4330cd1701db2fa249e45e27239b1a274ad24d`
