# Independent B5 dispersion implementation review

Reviewer: rf0_senior_review. Date: 2026-10-16.

Reviewed `ef0c77b63542944ae84c41bcfcf7935527d4d3ab` against `c84e85f`, plus the source-derived taxonomy review `/private/tmp/sceptre-dispersion-taxonomy-review.md`. **No P0/P1/material P2 estimator-state or test-transformation finding remains in this bounded implementation.** Final installed-source, supported-lane and workflow qualification remain open. No production edits, builds or tests were performed here.

## Production boundary

The only estimator expressions changed are `catch(...) {}` to two empty concrete catches for `std::domain_error` and `std::overflow_error`. `<stdexcept>` and an explanatory source contract are added. I verified the retained original is byte-identical to immutable Git21f9ea0, and independently compared comment/include/whitespace-normalized arithmetic/control/state text after reversing only that exact catch change; it matches.

The original whole-MLE/MM try boundary, pilot computation before it, estimate/method assignments, warning checks and result construction after it remain intact. A numerical MLE failure keeps pilot/method3; a numerical failure immediately before MM keeps the already committed MLE estimate/method1, including warning-flagged MLE state. No reset to pilot, new convergence policy or numerical fallback is introduced. Allocation, length, system, Rcpp and nonstandard failures now escape that boundary as intended. Pilot/result-allocation errors were already outside it.

The concrete catches agree with the inspected BH default binary64 CPU digamma/trigamma path: pole/domain and checked narrowing overflow are the reached throwing policies, fixed precision53 uses rational/asymptotic paths, and underflow/denorm default to ignore. MM has no Boost special-function call; ordinary floating NaN/Inf/errno outcomes are not silently converted into new C++ exceptions. This is not an assertion about every future BH version, custom policy/multiprecision configuration or possible exception taxonomy. Per-lane header/precision/policy evidence and original numerical parity remain required.

## Test transform and state proof

`test-dispersion-exceptions.R:20-33` fixes the original digest and rejects any candidate arithmetic/control/state difference other than the concrete catch substitution before adding instrumentation. The transform inserts one statement immediately before each unique actual MLE/MM invocation, retaining numerical functions and actual state assignments. Both original and candidate are compiled in separate namespaces in one test-only TU. The package has no hook, mutable injection API or substituted numerical implementation.

The probe's outer classifier is deliberately outside the actual estimator. Controlled domain/overflow injections compare exact original/candidate returns and reached stages; actual limit2 MLE executes its warning branch before the MM injection. The committed stage-reference crosscheck also compares all old blanket-catch results against the actual old numerical state at the same point. MLE pilot value and both method codes have independent assertions. Controlled allocation/system/length/Rcpp/unknown failures must propagate from candidate while old suppression/state is recorded. No real OOM or system mutation is used.

The 107 retained numerical inputs run through both actual bodies without injection. Exact result, propagation and stage equality are checked; this is stronger than comparing two handcrafted estimator copies. It is not exhaustive input coverage. The earlier log reports 107 numerical and 14 paired fault cases; exact committed runner rerun, including its added stage-state assertions, remains necessary before calling the committed proof passed. The runner accurately labels test-only compilation rather than installed-package qualification.

The platform fact probe records binary64 precision, evaluation type/digits and default policy values. CI runs it on Linux release/oldrel, macOS and Windows plus the dependency-free lane, and retains generated source/header/configuration evidence. Header hashes, package versions and compiler configuration are recorded, not interpreted as completed remote results. A lane with different dispatch/category must trigger source review rather than a broader catch or new scientific policy. No unresolved Eugene decision was identified.

## Reviewed hashes at commit

- `src/compute_nb_size_param_functions.cpp`: `27ea8665a984dcef83f35c646ab4b38989ea8cacafa28975aa7ffdedd7f9aeb2`
- `tools/refactor/dispersion-exceptions.cpp`: `f9c8a6e991adce64e09446c876ec4a25cf3bb5999328270de6453583124b29cc`
- `tools/refactor/test-dispersion-exceptions.R`: `6968a8278bba86da49529bdf1e0c63db4526ed6e29c5319be5d446b540a1fa6d`
- `.github/workflows/check-standard.yaml`: `039f48782766040c1e499eb13550c025e99ff130ee852afae86e084ae4499212`
- `tools/refactor/results/dispersion-characterization-21f9ea0/original.cpp`: `4add3e36820f9af6b794b715862c0981f44cb51f11750cf5107d39c2cd4ee776`
- `tools/refactor/results/dispersion-characterization-21f9ea0/report.rds`: `0d30f0ba5612e920a136e9f6236e2fef5469b6d4811c4b58de274db7835e6f2f`
