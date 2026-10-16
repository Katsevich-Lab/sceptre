# Independent CL2 allocation-remedy closure

Reviewer: rf0_senior_review. Date: 2026-10-16.

**The material P2 native-copy/R-longjmp finding is closed for the bounded count/mitochondrial kernels at `8841efa0a953fc545bf3b09197ab0e74e1968e1b`. No remaining P0/P1/material P2 found in this remedy versus ab477759.** This closes the specific amendment `/private/tmp/sceptre-cl2-allocation-review-amendment-ab477759.md` (SHA256 553cd8f20d74b78f1acc8f28e4a643937bcb16c8c6536e95d44c2683c78820e1), not generic package-wide Rcpp/OOM behavior or resource acceptance.

The six committed files matched immutable Git blobs. Retained successful probe copies of both actual native sources, generated bindings, probe and runner match this commit byte-for-byte. Untracked resource harness files are excluded. No source edits, compilation, failure injections, tests or benchmarks were executed by this reviewer; I inspected source, logs and retained basic RDS metadata.

## Ownership proof from source

- Generated bindings now convert only borrowed List/IntegerVector inputs for these three exports; no std::vector argument can be constructed before the body. Existing R-level formals/registration names and result schemas are unchanged, so R/RcppExports.R need not change.
- `sparse_counts.cpp:18-51,99-126,172-181` performs dimension/length/coordinate/sort validation, group coercion, input/ALTREP materialization and pointer capture, output/matrix allocation and named-list construction before native copies exist. Trivial pointer/length views use R_alloc. Newly normalized group lists and each coerced child are protected with Shield/PROTECT rather than PreciousPreserve, so no new normalized scratch is stranded through a skipped PreserveStorage destructor. The protected list retains each viewed R owner through computation.
- `sparse_counts.cpp:127-160,182-194` creates native coordinate/Boolean vectors only in explicit inner scopes. The helpers and loops use cached raw input/output pointers and pure C++ vector/arithmetic operations. No R API, Rcpp proxy, materialization, allocation/preservation, attributes or return-list work occurs inside. NT-sum overflow throws std::overflow_error so existing owners unwind before wrapper error construction. Native scopes end before return/copy/wrap.
- `sparse_matrix_functs.cpp:45-104` similarly creates all outputs, y scratch, full named result, materializes active mitochondrial/source inputs and caches all pointers before the sole owned mitochondrial vector. The old Rcpp-by-value loader call in that scope is replaced by the same ordered raw loads; ties, sums, loaded-value >0.5, repeated mitochondrial rows, disabled summaries and zero-total ratios are retained. Its native scope ends before return.

Coordinates still convert exactly once into owned zero-based vectors, retaining order/repetitions without caller mutation. NT overlap sums, complement totals, duplicate stored-coordinate presence, explicit zeros/fractions, empty shapes, pair ordering and integer capacity checks are unchanged. No ondisc path is changed. General Rcpp borrowed/output PreserveStorage behavior that already existed is not newly claimed universally leak-free; the new normalization buffers and native copies are the reviewed scope.

## Actual-source controlled evidence and independence

The probe interposes allocVector/allocMatrix/R_alloc/coercion/INTEGER/REAL/names/attributes/callable lookup before Rcpp template headers and again for the actual source includes. It observes precious preservation via the actual callable route, including return/copy/wrap preparation. Any R boundary encountered with a tracked native block is a failing safe C++ throw, not a skipped test or deliberate leak. A separate negative control deliberately visits such a boundary with a tracked owner and verifies rejection plus unwind to zero.

I independently checked that the retained instrumented source equals exactly (a) std::vector allocator replacement and (b) the single verified NT-total threshold substitution; no owner scope or calculation was moved by instrumentation. Nested vectors and vector<bool> allocations use the tracked allocator. The bounded threshold exercises the actual existing overflow throw without huge inputs. Each normal result is compared with a hand-specified independent count/covariate record and input serialization.

The completed pinned evidence has 214 main cases: 174 raw R-boundary errors, 18 bad_alloc injections, 18 native overflow injections, one bounded actual NT-overflow branch and three exact/nonmutation successes. All retained main-case states have live_blocks=0 and boundary_live_blocks=0. The separate unsafe-order negative control reports one rejected unsafe boundary and final live_blocks=0. Matrix allocation is observed through the installed Rcpp vector-plus-dim route; both that route and allocMatrix are interposed, without requiring an unexecuted route to pretend it ran. Earlier coverage-assertion-only failures remain separately preserved.

This is test-only actual-source evidence, not a new production injection interface, a real OOM experiment, or a universal claim for every Rcpp version. Header/normal package and supported-lane qualification remain the coordinator's gates.

## Remaining acceptance

The fresh immutable -O2 focused install/gates at `/private/tmp/sceptre-cl2-ownership-release-pp81irsy` were still ongoing during review. Complete them and qualify the final integrated immutable source with the normal frozen/native/storage/backed/cache/consumer/CI gates. Restart resource trials on this corrected pin; ab477759 resource results cannot qualify changed source. Triplet sorting workspace and common count resource limits remain unmeasured/open here; no correctness-cost exception is approved.

## Reviewed source SHA256

- `src/RcppExports.cpp`: `8057c466fd979687327623afb5813d48331a96a324764264c6178d0a71fc6909`
- `src/sparse_counts.cpp`: `60ac8650625c04e26b434e7b1447094bc9f40d3f584bc51d2e0e6744cff16711`
- `src/sparse_matrix_functs.cpp`: `a6368ee588f2f2c0b8f3fab54fc468bb5a1dfe93d972754a958944a5de255e99`
- `tests/testthat/test-refactor-storage-boundaries.R`: `34bd29c06e25d8128cdaafe37cfbf4e9a27c949f01377856883df50ffcd84ac6`
- `tools/refactor/count-ownership.cpp`: `e19111ea2ace8fba7490fb1aa6daec7797757a78d0f0152888b5ce1b8ea3c96f`
- `tools/refactor/test-count-ownership.R`: `1333ed0ed8b018f4afb5b5ea28216c99e8c6fba8667d5c71248fe794bdaea19b`

## Retained evidence SHA256

- `/private/tmp/sceptre-cl2-count-ownership-matrix-route/count-ownership-results.rds`: `2e219d54ae183293c215c374e3183a5e8ba7777386b7661459150aaf6f8b6b66`
- `/private/tmp/sceptre-cl2-count-ownership-matrix-route.log`: `cfebdeb4f2535846fb8a6b004161d7cc8bd752760de2d813a5370854832a1a7d`
