# Independent CL1 source review: ownership closure

Reviewed commit: **5ac0fca0063b5fc1ed7e0498de7658372d1fa526**. Focused revision comparison: dcba41858d6aa59cfb6eb056f8bc9c6d58eb2dc6 to the reviewed commit; original CL1 base remains 6a43fa5b2d32109a67e42138d3d845a6e5a2805a.

Outcome: **The prior material P2 is closed. No remaining P0, P1 or material P2 finding was identified in this focused source review.** This is source-review closure for CL1, not package/release readiness or a claim that pending immutable-build qualification has passed.

No source files were edited, tests/benchmarks rerun, or unsafe/OOM examples executed by this reviewer. The only write was this review artifact. Every file changed across the complete CL1 unit matches the exact reviewed committed blob. Frozen refactor fixtures have no diff from the CL1 base.

## Ownership finding closed

Each of the four current sampler producers now creates a protected NULL-address Rcpp external pointer and registers its original finalizer before its first native bank allocation. The fallible R external-pointer construction, precious preservation and finalizer registration therefore occur while there is no native bank to abandon.

During subsequent C++ bank construction, unique_ptr owns the vector-of-vectors, including its completed rows. A C++ allocation/other C++ exception destroys that bank through ordinary unwinding, while the protected external pointer remains NULL. The successful final transfer in sampling_safety.h uses only R_SetExternalPtrAddr and unique_ptr::release, followed by SEXP conversion; it does not allocate/preserve/register a new R pointer. At that point the existing finalizer owns the same bank type exactly once. No pointer class, bank layout, tag, persistence format or normal finalizer policy changes.

Inspection of the intervening sampler loops found no new fallible R API operation after native bank allocation and before transfer. Sampler arithmetic, loop/draw order and producer reuse scope remain as in the previously reviewed implementation. The specific late-registration R-nonlocal-jump leak is therefore removed, rather than merely converted to another exception path.

## Controlled R-error probe reviewed

The new developer-only probe includes the actual candidate generate_samples_functions.cpp and guard header, not a model sampler or frozen reference. Its translation-unit-local macro interposition reaches the three R boundary operations through Rcpp, and phase/error/address assertions are hard failures. It is run in a separate R process, introduces no production hook or installed-sampler replacement, does not request real OOM, and leaves the frozen native oracle untouched. The runner records pre/post input digests, generated instrumentation, dependency/session information and all case states, then retains source/probe copies.

The worker reports all 12 cases passed (four samplers times three boundaries). That runtime report was not independently rerun in this review. The probe demonstrates a NULL address at each injected R-error boundary; it is not a general live-native-allocation counter. Source inspection independently establishes that actual native bank allocation follows those boundaries. These two claims should remain distinct in qualification reports. The probe's avoidance of deliberately leaking an old non-NULL bank is appropriate.

## Other requested changes closed

- The public post-QC regression now omits calibration_group_size after reducing assignments to one active NT group while the unchanged metadata-derived default remains larger. It expects rejection and prevents a silent default clamp.
- R/analysis.R calculates and reuses one checked total for its permutation bank.
- The reached-stage skew-normal minimum is owned by the real fitter. The stage and internal quantile helper no longer repeat that guard. Stage-1-only short B2 metadata remains supported, while a reached short fit fails before its unsafe quantile access.

The previously reviewed valid-zero, checked-wide-arithmetic, post-QC authority and non-clamping semantics remain unchanged by this delta. Fixed-size inference still consumes only the required prefix of a sufficiently long row; hybrid suffixes are not validated against a smaller fit view; CRT consumes each actual possibly-empty row. All 17 frozen native inference cases and the existing independent full-bank gates remain in place. No per-pair full-bank scan or new statistical fallback is introduced.

## Remaining acceptance work

Fresh immutable-source installed-suite/native/original-oracle qualification is being produced separately for this exact commit. Preserve prior failed/stale-environment evidence and label regenerated same-environment oracle evidence by its actual source/environment. This review does not turn old test counts into current qualification. Subsequent cleanup units still need their own coherent migration/source/resource reviews.

## Complete CL1 changed-file SHA256 pins

- `R/analysis.R`: `775ed83ed0464cfb7bf754943d8c14db53b50d8fc0c8a6660139646f08d1a89c`
- `R/association_workers.R`: `50fcea5ebb2dd04de2e02f67dc71353a8882cb2b0a9c90ff966f708d0993adb2`
- `R/check_functions.R`: `331cf1755437f671439b18be5c9f49c3aaf9fb25818321c085f2518e9ec77246`
- `R/preprocessing_functions.R`: `3e4b4d94733325cf67da4c1263084bfc745c367d5c197a42558f60bb23a2ee4f`
- `R/s4_analysis_functs_1.R`: `326912e56118f33a443dac8eb07666b30fa2717060dd3e781b896f89de192629`
- `R/s4_analysis_functs_2.R`: `58ad5fff359ed6c12ededb41253ab0357c159d379c25766d18d7bd4a4c4113b7`
- `src/generate_samples_functions.cpp`: `23946d7a8a4721befef2f43c3670e33bc8d8b3485caaaaa23058319e10eb8642`
- `src/inference.h`: `71d284f711106435d62f9602872719476142f24d94a5ce3d880d29d86c533d66`
- `src/inference_scores.cpp`: `a05d3e221d37dcc05f83e4e403dc3a6ed7bc8a3536aa22fd67f877b458d3b0b0`
- `src/inference_stages.cpp`: `b447907d979ff90f1be75fd26fc1d93389973c8eb7e20c60f68d68bcdd16c5a4`
- `src/inference_tail.cpp`: `150551141185ea20e0e74f1b1c7595bf0683fa6c0b67f3b23e4bae17cb9206bd`
- `src/low_level_full_test.cpp`: `6a8819d33dc43266abd6bf86aaef1fec68cc4fdeea0761e88d70e003ecbd087d`
- `src/negative_control_functions.cpp`: `9a4a7d42fb667c9fc107849fb411166a3420a19af1b361fb1f3fc21bbc4eaf9f`
- `src/sampling_safety.h`: `0c801b9aa1ca7e90866109eab38851ae2c5e764300d3d1b02a8d9593afd879bd`
- `tests/testthat/test-refactor-calibration-safety.R`: `74378841f34b2736596e2cf64162f1a06870eb11121454a6126a71713987609a`
- `tools/refactor/sampling-ownership.cpp`: `a768902c1a7914a6ede8beb9f91b8e022ba0500ae4572708b4e1e4f5c71dfcea`
- `tools/refactor/test-sampling-ownership.R`: `c781899c36988a3e61657fb232c1a13a7eaae7fcd07f6b46847d46e2edbc1496`

Prior review retained: `/private/tmp/sceptre-cl1-senior-review-dcba418.md`, SHA256 `5d8ef72abee223ec7e743df88fce020e6111db2066b969c23edbb594f174d493`.
