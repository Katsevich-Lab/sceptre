# Independent pair/getter/mixture-contract reviews

Reviewer: rf0_senior_review. Date: 2026-10-16.

**No P0, P1 or material P2 found in the three bounded units below.** Reviewed immutable diffs and actual callers; no production edits, compilation, test reruns or benchmarks. Reported source-evaluated checks are not final installed/combined-source qualification. The separate symmetric-dense issue remains separately dispositioned and is not hidden by this outcome.

## B2/O1 — 2f40f4f

The cis kernel now evaluates the exact inclusive predicate across all supplied positions in supplied order, including repeated positions and later in-range blocks. int64 widening precedes subtraction/negation, so valid signed integer coordinates do not overflow; negative thresholds yield empty results. The NA/native-capacity guards establish representable native coordinates without inventing a biological missing-coordinate rule. Two passes allocate exact R output, avoiding an owned std::vector across R allocation/wrapping; local captures are trivial scalars. The new generated declaration matches IntegerVector with unchanged R formals/registration.

Midpoint sum widens before adding genomic endpoints; floor/integer midpoint semantics are retained when representable. Public unsorted/endpoint/wide-midpoint regressions use independent expected identities, and dangerous old arithmetic is never executed. Full traversal is necessary for unordered input, with O(n) work and only exact output allocation, rather than a full-length scratch or sort that changes order.

Trans exclusions copy the supplied positive-control table before setDT, so caller class/data/aliases remain unchanged. Product order, duplicates/exclusion rule, output class and expand.grid attributes are retained. Tests compare independent columns/class/nonmutation; the retained separate six-case full old/new comparison covers complete output attributes. The initial independent expected-table attribute mismatch is preserved and correctly excluded, not treated as a production failure/pass.

The retained candidate-source/native probe reports61 checks in3 tests, with source copies independently verified against2f40f4f. Its source environment resolves the compiled candidate predicate and candidate R constructors with helper auto-loading disabled. The code/help change stop→end matches the actual coordinate field and does not change API. Final immutable installed/native/pipeline qualification remains required.

## Assignment getter — 09a704c

The sole runtime increment_vector caller becomes unlist(...)-1L, creating the getter's owned zero-based slot while preserving valid integer source memberships/order, row identities, pointers and sparse class. Source/QC getter behavior, including the historical one-guide dimension drop, remains. Input/alias serialization and repeated calls are tested using valid increasing membership lists, with independent full/QC Boolean expectations, single/multiple/empty guide entries and exact slot checks. Reported24 source checks and6 exact old/new outputs remain targeted evidence; the invalid initial unsorted fixture is retained/excluded.

Current actual callers no longer need a native vector mutator. Updating its crosswalk to deletion rather than an unused renamed API is justified. Native export/registration retirement is a declared later coherent unit, not falsely completed here. Public formals/schema remain unchanged.

## Native mixture contracts — 5875e220

Independent normalized noncomment comparison confirms expressions unchanged. Prose agrees with the real R obtain_em_assignments caller: aligned count/Poisson-mean/log-factorial arrays, paired starts in order, no RNG/input mutation, retained arithmetic/label flip/stopping controls, strict greater-likelihood converged winner (ties retain earlier winner), zero-based winning start and initialized zero/FALSE/-Inf record when no winner. R owns active raw-count fallback. Relative likelihood tolerance retains the explicit -Inf case and ordinary other nonfinite/zero-denominator outcomes.

Minor precision suggestion: "mixture-density floor" is clearer than "probability floor", since code floors s=p0+p1 for log likelihood, not posterior Ti1s. This is prose only and no numerical regression run is warranted for unchanged expressions. Final source contracts/consumer qualification remain part of the coordinator's combined review.

## Reviewed immutable file hashes

### B2/O1 pair construction — `2f40f4f9b15262249a27aa710d408f70b4c6e04b`

- `R/pair_constructor_functs.R`: `95135d4871d73b90e53bba1f1f06cb3ae4210ddadc17155c2198ef7c685b7971`
- `man/construct_cis_pairs.Rd`: `745273cefdaf75e8076a48fa740ff4d5eeaa52a6131169daa84a89c6dd432193`
- `src/RcppExports.cpp`: `4582541709c7da2df03bd0e512ac054cc2ed7ed7b1e933b9d8068a3532b2b513`
- `src/pair_construction_functs.cpp`: `e6431b8a18ed986ee5f0a2c5e97ed05bc2a3729313960be00b974b76e8668038`
- `tests/testthat/test-refactor-pair-construction.R`: `6dab65733fd54a2a0ff8126894f0f760d5b96b496d04873693adad9de692ec70`
### Assignment getter — `09a704c968568866727cb061b1b8a9a83d6d13db`

- `R/getters_and_setters.R`: `bbbbbe812200e74218149c7b21f9ddf21bbfe915d371cde5458c3f04cb153e26`
- `tests/testthat/test-refactor-assignment-getter.R`: `7313c229b749267e9a724d81406f805d8f925f99242945e024c5cdbb52f69d7c`
- `tools/refactor/cleanup/symbol-dispositions.md`: `2a87d6990a86d8395e05bd51c7e60dd72d599029a654ce23fcab9e40986ae8d8`
### Native mixture contracts — `5875e220cfcfc70e444b6455abf0c645a49297ff`

- `src/mixture_functs.cpp`: `bd1051e6b353279c5f6ba85d87b8baec8c85f7b385545ebb7da39967fc10edb4`
