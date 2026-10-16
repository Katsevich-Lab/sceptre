# Independent B5/T12 dispersion exception taxonomy review

Outcome: **The pinned default CPU/double implementation supports replacing catch(...) with two concrete empty catches: std::domain_error and std::overflow_error.** Preserve the existing try boundary and estimate/method state exactly. No unresolved scientific-policy conflict was identified for that bounded replacement. This is a fact/recommendation review before production changes; controlled fault-state/propagation proof, final-source review and supported-lane qualification remain required.

No estimator tests, builds or production edits were performed. A lightweight R metadata query located the installed BH headers. The original source, transformation, report/probe/runner, default policies, precision dispatch and relevant header call chain were inspected.

## Evidence verified and its limits

- original.cpp is byte-identical to 21f9ea098b69c444a884a49bc254de711968e3b5:src/compute_nb_size_param_functions.cpp.
- All ten files covered by artifact-manifest.json match their SHA256 entries. The three header digests in preparation.json match the installed headers read here.
- The original-to-instrumented diff only removes the export, renames the entry and records the caught exception. Arithmetic, loops, fallback decisions and returned estimate/method are unchanged. The extra macro name in probe.cpp does not introduce candidate arithmetic.
- The retained report contains 107 aligned cases with original_exact=true: 61 domain errors, 35 overflow errors and 11 no exception. Every recorded numerical throw retains method3/pilot state. Those are reported executed facts of the retained run, not a rerun by this reviewer or an exhaustive input taxonomy.
- The independent installed original is selected explicitly by run.R and its actual package path is recorded. Final implementation evidence should also bind the reused immutable original build/DLL manifest and before/after hashes, rather than treat a library directory string alone as binary provenance.

## Source-derived numerical categories

The estimator calls digamma/trigamma with ordinary double arguments from th and th+y[i]. For the installed BH1.90.0.1 headers, default policies::precision<double,Policy> derives the tag from double's precision, not from the promoted evaluation type. A supported binary64 CPU lane therefore selects tag53 even if calculation is promoted to long double. Promotion and platform long-double range can change which inputs overflow, but not this fixed-precision dispatch family.

| Reached mechanism | Default behavior and typed catch |
| --- | --- |
| Digamma/trigamma poles/reflection | raise_pole_error forwards to the throw-on-error domain policy, which throws std::domain_error. Boost wrapping remains derived from that concrete type. Preserve suppression with catch(const std::domain_error&). |
| Checked narrowing to double | Overflow policy defaults to throw_on_error and raises std::overflow_error. Preserve suppression with catch(const std::overflow_error&). |
| Narrowing underflow/denormal results | Default policies are ignore_error and do not throw. Retain the existing zero/denormal/nonfinite arithmetic; do not add a new underflow fallback. |
| Fixed tag53 rational/asymptotic calculations | No series/root convergence exception path is selected. Trigamma's tag0 polygamma and digamma's generic series branches are not used for these double calls. |
| Reflection helpers/constants | sin_pi uses scalar arithmetic and ignores narrowing overflow; pow<2> uses multiplication. Built-in floating constants use direct numeric conversion, not lexical conversion. No additional exception category is required by this call chain. |

BH declares evaluation_error and rounding_error as generic policy categories; their declaration is not evidence that this estimator reaches them. Do not add those catches solely because they exist or were included in the diagnostic classifier. Likewise no generic std::exception, std::runtime_error or catch(...) is needed. Allocation while formatting/wrapping a numerical exception may itself fail with bad_alloc; that unrelated failure should propagate, not be reclassified as the intended numerical exception.

This conclusion is bounded to the inspected BH/default-policy CPU dispatch. It is not a claim for every historical/future BH implementation, custom global BOOST_MATH_* policy, GPU/multiprecision path, or disabled-exception configuration.

## MM and exact state preservation

MM calls ordinary scalar arithmetic, Rcpp sugar sum/pow and constructs its returned std::vector. Installed Rcpp sugar pow calls ::pow; sum performs the existing floating additions. NaN/Inf, division by zero and convergence flags are arithmetic outcomes, not evidence of a new C++ numerical exception category. MM has no digamma/trigamma call. Do not swallow its Rcpp, allocation, length or system failures under the label numerical.

Keep both typed catches at the original whole-MLE/MM try boundary, with empty bodies. In particular:

1. A numerical throw before MLE returns leaves the pilot estimate/method3, as all96 observed throws do.
2. If a numerical throw were to occur while evaluating MM before its result assignment, the original catch preserves the previously committed MLE estimate/method1, even if that MLE raised its numerical warning flag. Do not reset to pilot, rename the method or normalize state. This partial-state rule needs controlled original-versus-candidate evidence; the107 probes do not establish it experimentally.
3. Successful MM and its existing warning-driven pilot fallback retain their present assignment order and values.
4. Pilot computation and final List::create are outside the original try; leave those boundaries unchanged.

std::bad_alloc, std::length_error, std::system_error, Rcpp exceptions and unknown/nonstandard exceptions then propagate naturally. Ordinary floating-point status/errno behavior stays unchanged. This selects no estimator or convergence policy.

## Required implementation/qualification evidence

- Use test-only interposition of the actual original and candidate MLE/MM call boundaries, with source/transform hashes. Inject domain/overflow at MLE and MM points and compare exact returned estimate/method state; execute a real warning-driven MLE-to-MM path rather than relabel an MLE-only case.
- Inject unrelated bad_alloc and system/length/Rcpp failures at both boundaries without real OOM or production hooks. They must escape the candidate; the comparison should explicitly record the intended difference from the old blanket catch. Preserve actual numerical and ordinary cases independently of injected tests.
- On Linux release/oldrel, macOS and Windows/Rtools, record BH/header hashes, compiler flags, numeric_limits<double>::radix/digits, effective precision tag, promoted evaluation type/digits, and default pole/domain/overflow/underflow/denorm policy values. Normal GCC/Clang/Rtools binary64 configurations follow the tag53 argument above, but this review did not run those remote lanes. A test-only policy-fact probe can confirm the assumptions without changing production policies or adding a fallback catch.
- If any lane has a different reachable numerical category or changed dependency dispatch, review that concrete path and preserve its original result/state before adoption. Do not infer it from the 107 counts or broaden to std::exception to make a lane pass. No production policy assertion/minimum-dependency change is implied by this recommendation.
- Bind final-source immutable candidate/original binaries and rerun the exact numerical cases, normal scientific fixtures and relevant platform tests. Passing taxonomy probes alone does not qualify the replacement.

The current remaining work is evidence and implementation review, not an unresolved user choice. If actual characterization exposes a new numerical fallback/state decision, stop that dependent correction with its precise conflict as the approved cleanup plan requires.

## Evidence artifact SHA256

- `tools/refactor/results/dispersion-characterization-21f9ea0/artifact-manifest.json`: `5b4c93ac26ef077c267d2d596f66413bf7e42f6d94d70ece8f04806f7ff23fb7`
- `tools/refactor/results/dispersion-characterization-21f9ea0/failed-initial-loader.log`: `0a06b6f1d34de86bdb8b2b86a107fcd47537ea013c063b3700257a8696ee899e`
- `tools/refactor/results/dispersion-characterization-21f9ea0/instrumented.inc`: `7321a095bec46a3111169cf1a6fa9363aa1355d5b81726db1692272434db7565`
- `tools/refactor/results/dispersion-characterization-21f9ea0/original.cpp`: `4add3e36820f9af6b794b715862c0981f44cb51f11750cf5107d39c2cd4ee776`
- `tools/refactor/results/dispersion-characterization-21f9ea0/preparation.json`: `f459bd6571e499b11a8ffec7f68b6e51fa915869affb515ba2c7cc977b66c4d8`
- `tools/refactor/results/dispersion-characterization-21f9ea0/probe.cpp`: `c337b3707188533483373714cdcf551a854c810d080083c2c4a391506b76bd2d`
- `tools/refactor/results/dispersion-characterization-21f9ea0/report.json`: `127f4d5b72ec94da1dfe700449615b9f8f60bad98802ce738ebeab33203587f6`
- `tools/refactor/results/dispersion-characterization-21f9ea0/report.md`: `d78c67c37a8d507acdb11fdec3025fbac82d762356f9f558e900c2bf67490d99`
- `tools/refactor/results/dispersion-characterization-21f9ea0/report.rds`: `0d30f0ba5612e920a136e9f6236e2fef5469b6d4811c4b58de274db7835e6f2f`
- `tools/refactor/results/dispersion-characterization-21f9ea0/run.R`: `d54ef4840f07d97c0f23c0f83d285f7069c69bb9a9a6c512c1d72b84b7600b02`
- `tools/refactor/results/dispersion-characterization-21f9ea0/run.log`: `c5fadfd49748aa038949ff979842995d3d2320f84ed8b5baf72416c2ca4a51d1`

## Inspected installed header SHA256

- `/Library/Frameworks/R.framework/Versions/4.6/Resources/library/BH/include/boost/math/policies/policy.hpp`: `0f0006710174f5bd9bacb6bdd860cf512282124f19c621877741822c9c135375`
- `/Library/Frameworks/R.framework/Versions/4.6/Resources/library/BH/include/boost/math/policies/error_handling.hpp`: `5e814d688c6450caa45e4e18934dd9e713d49204800be7829b8f9d5091096e03`
- `/Library/Frameworks/R.framework/Versions/4.6/Resources/library/BH/include/boost/math/special_functions/digamma.hpp`: `9c51bf9e313062f50cb9b21e2a29e9c9210f27eb9689003b38a40e2552f9cb0f`
- `/Library/Frameworks/R.framework/Versions/4.6/Resources/library/BH/include/boost/math/special_functions/trigamma.hpp`: `b1ae9e027ad15523504cc669a2f6de18edf0920cb3e3be9276923ed4c64d1ad4`
- `/Library/Frameworks/R.framework/Versions/4.6/Resources/library/BH/include/boost/math/special_functions/sin_pi.hpp`: `b0b6ef3085e697a36f68e61e94dcbeaf511b74c285c16a6a5e4034205a851d31`
- `/Library/Frameworks/R.framework/Versions/4.6/Resources/library/BH/include/boost/math/special_functions/pow.hpp`: `498a9052d42938fb2e239a1520ff71bbfa45cd7479cc0127170fdb8a762ef820`
- `/Library/Frameworks/R.framework/Versions/4.6/Resources/library/BH/include/boost/math/tools/big_constant.hpp`: `cc1005b66eaf094d251b2c7185e7aa118d546c7ccb1cb61416725ab5d996de84`
- `/Library/Frameworks/R.framework/Versions/4.6/Resources/library/BH/include/boost/math/tools/throw_exception.hpp`: `6f0ebc04d8c47e30f48718bd8ffdae4f2931ab23c9371b275832a950b7a07409`
- `/Library/Frameworks/R.framework/Versions/4.6/Resources/library/BH/include/boost/math/tools/numeric_limits.hpp`: `bda01f379ea4fad40676bf4c4d68eaf5caa884833df38ad6cd2f619f1f4a7793`
- `/Library/Frameworks/R.framework/Versions/4.6/Resources/library/BH/include/boost/math/tools/config.hpp`: `937bfeba8cc4655b6e77b458d8d7fd0efca1a5dc51463f02a60406ffc1d6639d`
- `/Users/ekatsevi/Library/R/arm64/4.6/library/Rcpp/include/Rcpp/sugar/functions/pow.h`: `54c337000aaf58cac71b539abaa525e41e5659a16fd2e27698f2336603ad851d`
- `/Users/ekatsevi/Library/R/arm64/4.6/library/Rcpp/include/Rcpp/sugar/functions/sum.h`: `f07b6efe76aef1e924b6ee77c8e84e2dbbf3adf12770484fdf3173c86287cbed`
