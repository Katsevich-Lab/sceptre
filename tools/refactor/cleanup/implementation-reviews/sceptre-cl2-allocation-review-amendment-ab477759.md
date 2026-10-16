# CL2 amendment: R longjmp bypasses newly owned native vectors

Reviewer: rf0_senior_review. Date: 2026-10-16.

**Material P2 — corrective work required before CL2 acceptance/resources.** This amendment supersedes the allocation-safety aspect of `/private/tmp/sceptre-cl2-senior-review-ab477759.md`. Its coordinate/count/ordinary ownership assessment remains valid. The first pass did not inspect installed Rcpp's raw-R allocation unwind behavior; that missing inspection matters and is corrected here. Root accepted this finding and held resource trials. No production edit, failure injection, test or benchmark was executed by this reviewer.

## Source-derived finding

Affected source is `ab4777591cb31906f1c154c1825f0a6bed8c7c6a`:

- `src/sparse_counts.cpp:65-88,115-116`: exported nested group/request/cell `std::vector` parameters are owned when Rcpp IntegerVector/IntegerMatrix outputs allocate. The named return List/names are also built while both coordinate owners and Boolean scratch remain alive.
- `src/sparse_counts.cpp:125-133`: exported group/cell owners precede IntegerMatrix allocation; a return-copy/preservation operation must also be outside native-owner lifetime.
- `src/sparse_matrix_functs.cpp:44-48,93-95`: the new owned mitochondrial vector precedes output/y allocation and named return-list construction.
- `src/RcppExports.cpp:203-233,264-278`: generated parameter adapters eventually call `InputParameter<T>::operator T()`, constructing owned vector arguments for the function call before its body. Merely moving allocations to the start of the current body leaves these parameter owners alive.

Installed Rcpp Version: 1.1.1-1.1 directly uses `Rf_allocVector` in Vector constructors and `Rf_allocMatrix`/dimension allocation in Matrix construction. PreserveStorage also invokes PreciousPreserve. Its BEGIN_RCPP/END_RCPP wrapper is C++ try/catch; there is no automatic R_UnwindProtect around these allocations. A raw R allocation error can longjmp back to R, skipping C++ vector destructors. The resulting native malloc storage is not reclaimed by R GC. Count group/cell copies can be large, so recoverable output-allocation failures can leave substantial permanent process memory until exit. This is a new native-copy regression even though old Rcpp and Boolean-scratch lifetimes already had related generic failure-path limitations. No ordinary count/scientific result regression is asserted.

The optional `Rcpp::unwindProtect` helper does not justify an easy blanket wrapper: it allocates a continuation token before its setjmp, and an R longjmp from a callback skips callback-local C++ owners before it is rethrown as LongjumpException outside that callback. Only already-existing caller-side owners are then unwound. A whole-compute callback therefore does not repair inner owner lifetimes.

## Bounded required remedy agreed with native worker

1. Export borrowed R inputs (List/IntegerVector or SEXP), eliminating generated owned-std::vector argument conversion before the function starts.
2. Before any owned native coordinate/Boolean vectors exist: validate dimensions/lengths/coordinates/sort, perform required coercion, materialize ALTREP data, capture pointers/lengths, allocate R output vectors/matrices/dimensions and fully named return lists, and finish R protection/preservation work. Keep normalized owners protected for the whole computation. Trivial pointer/length views in checked R_alloc storage are appropriate and avoid repeated potentially allocating data-pointer calls while native owners are alive.
3. Create owned native copies once from those cached pointers and compute in a bounded inner scope using raw input/output pointers only. No R API allocation, Rcpp proxy coercion, R callbacks, output/list/name creation, attribute assignment, new data materialization or preservation may occur there. Native allocation failures and NT-sum overflow must throw C++ exceptions that unwind owners before wrapper error conversion; no raw R error in that scope.
4. End the native scope before return/copy/wrap/preservation. Preserve existing indices, result dimensions/field names/order, count multiplicities, stored-zero semantics, native nonmutation and retained ondisc copies. No production hook or generic allocation framework is required.

Protect NEW normalization scratch with R PROTECT/Shield (or NoProtectStorage plus Shield), rather than adding PreciousPreserve-owned copies whose destructor could itself be skipped on a later preflight R error. Existing generic Rcpp wrapper/precious-object behavior should not be mislabeled universally solved by this bounded correction. The acceptance claim should describe the new native-copy/scratch lifetime guarantee precisely.

## Required controlled proof

Use test-only interposition of actual candidate source, with transformation/source hashes and no real OOM. Trigger raw R errors at R vector/matrix allocation, normalization/coercion/materialization, R_alloc/view preparation, List/name/attribute/preservation phases as applicable. R catches must observe zero live tracked new native-copy/scratch allocations. Exercise actual compute allocation/overflow exceptions separately and prove owners unwind before error reporting. Normal calls must retain existing exact count/shape/nonmutation tests and frozen scientific parity. Test-only allocation/owner accounting must correspond to actual candidate owners; do not prove a stand-in implementation. Rebuild and qualify the resulting new immutable source, then restart resource trials on that pin; ab477759 timings cannot qualify changed code.

## Installed headers inspected

- `/Users/ekatsevi/Library/R/arm64/4.6/library/Rcpp/include/Rcpp/InputParameter.h`: `8a2894ce2f8a43547c24b7ff92f74272c3e103f8a8fe6b1b0082a9a55e4d74d0`
- `/Users/ekatsevi/Library/R/arm64/4.6/library/Rcpp/include/Rcpp/vector/Vector.h`: `9da784c28970aefee8e0f8904486b4f2b9ce2e1bff48d08047ab7c9aef5011ca`
- `/Users/ekatsevi/Library/R/arm64/4.6/library/Rcpp/include/Rcpp/vector/Matrix.h`: `b9f48c70a8cdfd8370b6879b90fb05d24f99302f4da8fb6fee1a393e7b130ce5`
- `/Users/ekatsevi/Library/R/arm64/4.6/library/Rcpp/include/Rcpp/storage/PreserveStorage.h`: `147e22721e9bfd810a3a65f055489aaa69c816f90815b4ed9efdae5b906fa238`
- `/Users/ekatsevi/Library/R/arm64/4.6/library/Rcpp/include/Rcpp/macros/macros.h`: `93c4d6e0520ee3975f970e726cc0685eee9597d80ef7b7e3372aa67d34ee4b9a`
- `/Users/ekatsevi/Library/R/arm64/4.6/library/Rcpp/include/Rcpp/unwindProtect.h`: `f278400619765a48c99adade90cf448cb3f7f2da9daac13c0b0884924c728fd1`
