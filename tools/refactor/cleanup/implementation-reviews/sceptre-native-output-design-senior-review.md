# Independent native output lifetime design review

Reviewed base: dc093d27acefd2b3897eaf8b5a0bb4e9cbe7282f. Reviewer: independent senior source reviewer. This reviews the proposed exact-result conversion design, not an implementation or runtime qualification. No R/compiler/test process or production edit was performed for this review.

## Outcome

Approve the bounded design to implement. No remaining material design defect identified in the preprovided-token transfer sequence. The existing native lifetime findings remain open until actual source and controlled probes are reviewed.

The group allocation-first path removes R allocation from the lifetime of its set/sample owners without eager sizing. Pair and inference outputs retain exactly the accepted/selected native result until a protected, exact-length conversion; neither rejected-pair capacity nor an unreached stage budget is allocated.

## Verified boundary and required implementation invariants

R4.1 context.c lines914–955 shows that a preprovided non-NULL continuation bypasses the allocation branch. begincontext (lines247–280) contains context/stack assignments and vmaxget, followed by SETJMP before entering the leaf. The cleanup callback executes after endcontext and before R_ContinueUnwind. On success the result is assigned to CAR of the protected continuation before cleanup. These facts support the proposed transfer without adding a C++-exception bridge across R frames. The same ordering is present in the recorded current R source. Public R4.1 Rinternals.h declares both APIs; DESCRIPTION already requires R>=4.1, Rcpp>=1.0.9.

1. Protect the continuation before any native result owner. Balance that protection on ordinary and C++ compute-error paths.
2. Cache/materialize R inputs before native computation. End the complete compute scope normally, including native scratch, unique_ptr state, observed-position Rcpp owner and XPtr wrapper. Transfer only the raw payload pointer and trivial fields; no R allocation, preservation, materialization or conversion may intervene before R_UnwindProtect.
3. The protected callback is an exact-size raw-C-API conversion leaf. It has no C++/Rcpp precious owners, allocating C++ operation, C++ throw, user callback, installed handler, or Rcpp::stop. Raw R allocation/errors are handled by the registered cleanup.
4. Cleanup is noexcept and allocation-free, deletes the payload exactly once on normal and R-jump exit, and leaves no payload for the generated final wrapper.
5. Preserve result names, scalar types, distribution order/cardinality, selected stage semantics, RNG and loop order. Keep the result rooted via the protected continuation through cleanup until raw SEXP return.
6. Probe actual source scopes, not a rewritten surrogate; inject at observed token/input/output/name/attribute/return boundaries and native allocations, preserve failures and safe negative control, and assert payload/scratch release without a real OOM. Check ordinary outputs independently.

## Qualified scope

This approves native payload/scratch ownership and an owner-free conversion leaf. It does not certify universal Rcpp borrowed-parameter precious-node destruction, generated-wrapper cleanup, every package native lifetime, interrupt safety, or arbitrary OOM behavior. A raw R jump can still bypass unrelated generated C++ parameter wrappers; that broader existing boundary is not silently blessed by this repair. Actual source must avoid adding callback-local precious scratch. Runtime/component/frozen/installed gates remain required and separate.

## Primary source references

- R4.1 context implementation: https://raw.githubusercontent.com/wch/r-source/R-4-1-branch/src/main/context.c
- R4.1 public declarations: https://raw.githubusercontent.com/wch/r-source/R-4-1-branch/src/include/Rinternals.h
- Current R context implementation: https://raw.githubusercontent.com/wch/r-source/trunk/src/main/context.c

## Reviewed byte pins

- `/private/tmp/sceptre-native-output-lifetime-design.md`: `ff04621e049d81c219f267c628447d1d37c2281800167b7d92201fae7e379466`
- `/private/tmp/sceptre-native-output-design-sources/R-4.1-context.c`: `f76ddb376b6c99896d97bfa9901c2d3fc739ce5e4c4eb64ac3caefcd45a39015`
- `/private/tmp/sceptre-native-output-design-sources/R-4.1-Rinternals.h`: `f1568a37814dbc1ff6d687dce368cbbbb73029dc85fa8a85c06213f3cb3c2ab1`
- `/private/tmp/sceptre-native-output-design-sources/R-current-context.c`: `4b0c535df65d8e0caefc7e26e52d20f630663969f2cf60035f4b0e8b4bc534d5`
