# Independent CL1 implementation review

Source reviewed: dcba41858d6aa59cfb6eb056f8bc9c6d58eb2dc6 against 6a43fa5b2d32109a67e42138d3d845a6e5a2805a in the calibration-safety worktree. Outcome: **CHANGES REQUESTED: one material P2 ownership gap.** No P0/P1 implementation finding was identified in this focused review. This is source/design review; test counts are not used as proof of safety or readiness.

No production files were edited, tests/benchmarks run, or unsafe examples executed. A lightweight R metadata query located the installed Rcpp headers, which were read to resolve the registration exception behavior. All 15 changed files match the reviewed committed blobs; the frozen fixture directory has no diff in this commit.

## Material finding

### P2 — R errors during late pointer registration still bypass bank cleanup

Location: src/sampling_safety.h:46-49, called after bank construction by src/generate_samples_functions.cpp.

The unique_ptr now correctly destroys partially constructed banks when C++ allocation or other C++ exceptions unwind. However register_bank first allocates/preserves an R external pointer and registers its finalizer only after the entire bank exists. Installed Rcpp XPtr.h calls R_MakeExternalPtr directly, then PreserveStorage::set__ calls Rcpp_PreciousPreserve, and setDeleteFinalizer calls R_RegisterCFinalizerEx. These are not enclosed in unwindProtect by XPtr or the generated BEGIN_RCPP wrapper. An R allocation error/nonlocal jump in this registration sequence bypasses the C++ unique_ptr destructor, leaving the full bank allocated without a registered owner. Therefore the planned S3 R-error path is not fully closed; normal success and std::bad_alloc tests cannot establish it.

Concrete correction: create/protect a NULL-address external pointer and register its normal finalizer before allocating/building the bank. Keep the bank in unique_ptr during construction; once complete, install its address and release unique ownership with nonallocating operations. Alternatively explicitly use Rcpp's supported unwind protection around the fallible registration path and verify exact transfer/no double ownership. Keep the same vector-of-vectors type, pointer representation and finalizer. Review the failure path and, where practical, use a controlled test-only R error rather than real OOM to exercise the chosen ownership transfer. No new pointer class/schema is required.

## Verified source properties

- Public calibration validates group size against actual post-QC nonempty NT groups before combinations. It retains the old metadata-derived default and rejects an oversized default/request; it does not silently change the calibration treatment definition.
- Complement calibration clamps only its conservative sum upper bound to the cell count. The formerly safe M values, sampler bodies, generator families/seeds, loop order and returned bank order are unchanged.
- Fisher-Yates safe zero domains and CRT zero-size/endpoint domains are retained. Hybrid requires positive N/m, checks N+M and M+1 before their int expressions, and retains the same prefix-producing loop.
- The fixed-size native consumer requires row length >= n_trt and reads only that prefix. CRT uses each row's actual length, including zero. No full-row equality rule or smaller-fit-view scan of an unused hybrid suffix is introduced. No full bank-content scan runs per pair.
- Score dimensions and observed count/ranges are checked before observed indexing. Budgets and selected segment ends use checked wide sums. The reached-stage SN check leaves an unused short B2 alone, and B3=0 fallback segment/record behavior remains unchanged.
- Candidate-pair products and treatment sums use wide arithmetic before narrowing; nonnegative input counts plus treatment<=total make the final int subtraction safe. Group enumeration/sampling check actual bounds, represented combination count and output matrix capacity. Finite possible-group bounds clamp the double estimate before narrowing; +Inf remains only an upper bound, never a cast/allocation count.
- R-derived B3 uses the existing double multiplier/ceiling arithmetic and checks the complete sum before narrowing. No new rate floor, streaming sampler, scientific fallback or arbitrary cardinality cap was added.
- All 17 existing native inference cases and frozen input/source/expected artifacts remain unchanged. Exact committed-build runtime qualification is being produced separately and was not rerun here.

## Nonblocking cleanup and coverage improvements

1. R/analysis.R:14 validates/discards checked_resampling_total and calls it again in the permutation branch. Compute once and reuse the result; the later context migration can carry one authoritative checked budget.
2. The >=100 SN invariant is checked in evaluate_staged_tail, fit_and_evaluate_skew_normal and check_sn_tail. Prefer one real owning fitter boundary, or explain any retained caller-specific guard, rather than leave repeated ownership after the planned native cleanup.
3. The public post-QC test covers an explicitly oversized group size. Add the corresponding omitted/default group-size regression when metadata-derived default exceeds the surviving group count, confirming no silent clamp. Existing frozen native tests retain all hybrid prefixes and variable CRT cardinality; keep those gates active.

These improvements do not introduce a second material blocker. Final safety/resource/source approval depends on the corrected immutable source and its required qualification; this review does not infer readiness from 341 focused checks.

## Reviewed changed-file SHA256

- `R/analysis.R`: `dcd3f754fb7840afbbd8dbf353ede0fb8a18e184a45b857830bc6447416dbbbb`
- `R/association_workers.R`: `50fcea5ebb2dd04de2e02f67dc71353a8882cb2b0a9c90ff966f708d0993adb2`
- `R/check_functions.R`: `331cf1755437f671439b18be5c9f49c3aaf9fb25818321c085f2518e9ec77246`
- `R/preprocessing_functions.R`: `3e4b4d94733325cf67da4c1263084bfc745c367d5c197a42558f60bb23a2ee4f`
- `R/s4_analysis_functs_1.R`: `326912e56118f33a443dac8eb07666b30fa2717060dd3e781b896f89de192629`
- `R/s4_analysis_functs_2.R`: `58ad5fff359ed6c12ededb41253ab0357c159d379c25766d18d7bd4a4c4113b7`
- `src/generate_samples_functions.cpp`: `993cf24473dd21f7f9403dd5c24e7792aa61fc8fe95cfc2343edfc85df594cf8`
- `src/inference.h`: `71d284f711106435d62f9602872719476142f24d94a5ce3d880d29d86c533d66`
- `src/inference_scores.cpp`: `a05d3e221d37dcc05f83e4e403dc3a6ed7bc8a3536aa22fd67f877b458d3b0b0`
- `src/inference_stages.cpp`: `725c25a559758cbf846281e008e0fabf8cfdb2857417f86fdb39a66050b952db`
- `src/inference_tail.cpp`: `cca4f3c4f35164c6142d7a7d74205f43e96e79758624dfd8fd277d9fb3184469`
- `src/low_level_full_test.cpp`: `6a8819d33dc43266abd6bf86aaef1fec68cc4fdeea0761e88d70e003ecbd087d`
- `src/negative_control_functions.cpp`: `9a4a7d42fb667c9fc107849fb411166a3420a19af1b361fb1f3fc21bbc4eaf9f`
- `src/sampling_safety.h`: `bc02aab6c030d426fdc3bafb7e90761a5d23332dcbe23b7fecd40b32a4290c3b`
- `tests/testthat/test-refactor-calibration-safety.R`: `68e61bfc676270669bd49c66446f1115ea7d16147fa0edf5e6aac5dcc11a2749`

## Rcpp source evidence

- `/Users/ekatsevi/Library/R/arm64/4.6/library/Rcpp/include/Rcpp/XPtr.h`: `2c9dd5b0ff4be4fd9e4bd9d87a683877216d078b10c5ffddd646232dca7700f8`
- `/Users/ekatsevi/Library/R/arm64/4.6/library/Rcpp/include/Rcpp/storage/PreserveStorage.h`: `147e22721e9bfd810a3a65f055489aaa69c816f90815b4ed9efdae5b906fa238`
- `/Users/ekatsevi/Library/R/arm64/4.6/library/Rcpp/include/Rcpp/macros/macros.h`: `93c4d6e0520ee3975f970e726cc0685eee9597d80ef7b7e3372aa67d34ee4b9a`
