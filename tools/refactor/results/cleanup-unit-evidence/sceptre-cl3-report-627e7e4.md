# CL3 implementation — 627e7e4

The coherent source commit is `627e7e4c619b06359980b9e04fb246faa67c398d`, after the count-lane commit adopted as `f8f55b1`. The latter is the cherry-pick of `ab4777591cb31906f1c154c1825f0a6bed8c7c6a`; cherry-pick only CL3 when integrating onto a parent that already has that count unit. The JSON companion pins changed source and DLL bytes, build paths, logs and settings.

## Responsibility changes

- `analysis_context.R`: ten documented named functions own coordinate maps, saved NT encoding, request resolution, one lean treatment closure and one common matrix design/view. NT encoding uses actual membership-to-pool matches; valid disjoint saved memberships/order remain. Context has required eligible pairs and checks the full stage budget before preparation. Saved covariate design remains full-source, and fit-cache strict full-matrix authority is documented accurately. No per-pair full views or all-guide view/design bank.
- `analysis.R` and the three facade callers: direct explicit permutation/CRT dispatch through the prepared context and ordinary execution settings. One complete permutation bank is generated before traversal. Full-family result assembly remains. Public formals and S4 slots are unchanged.
- `association_execution.R`: separate response-oriented permutation and guide-oriented CRT drivers, with three documented block workers. Prepared integer response/group coordinates replace repeated ID rematching. Existing first-request/seeded block order, requested row order, shared fit insertion, warm-cache row reads and CRT prepass then guide-phase append lifecycle remain. One cache lookup per existing fit hit. NT guide fits remain uncached.
- `association_workers.R`: four explicit required-input cases retain scientific scope. Shared permutation score/fit reused over groups; NT permutation fits per group. CRT fits assignments and draws one full B1+B2+B3 bank per guide before its response loop; shared response fits come from the cache, NT response fits remain per guide. Existing native names remain for later coordinated retirement.
- `execution.R`: six documented helpers keep dispatch, table binding, file lifecycle, reporting and maintenance distinct. ID lookup is visible in the caller. Progress remains item 1/every 5, with exact serial/parallel newline bytes; GC remains every 200 before read/fit and independent of display.
- `fitting.R` / `score_preparation.R`: direct names `fit_response`, `fit_assignment_probabilities`, `prepare_score_components`; unchanged full-statistic formulas and response fitting. Uncalled coefficient-only and reduced-score modes and their three forwarding aliases are deleted, including the empty precomputation module.
- `mixture_model_functs.R`: documented guide-block worker uses explicit progress/GC. Remove the uncalled, broken false initialization mode; retain the existing Poisson initialization, warning handling, native EM success and real low-count/nonconvergence fallback. Fixed starting guesses/caller RNG, order, diagnostics, dispatch and log bytes remain.

## Named behavior correction and evidence

B3 alone supersedes the old NT one-column matrix-drop error. Matrix dimensions are preserved at the common/group boundaries. Four independent intercept cases cover calibration/discovery and permutation/CRT, with source positions written independently of context helpers and an independent scalar intercept-score calculation. The historical edges.rds artifact and original-source edge-capture helper remain unchanged; only the live intercept_only_nt error assertion is superseded.

The new tests additionally assert seeded factor-label/integer partition equivalence, exact progress/file bytes, fixed maintenance cadence, public cold/warm fit/bank/row-read scopes, actual native EM success and nonconvergence, low-count fallback, caller RNG, serial/fork mixture assignments/order and two 201-item progress logs. Existing native/index/cache/fitter/saved/API/backed/serial-fork/scientific fixtures remain active.

## Executed validation

Immutable candidate archive/install at `/private/tmp/sceptre-cl3-source-627e7e4` and `/private/tmp/sceptre-cl3-library-627e7e4`; clean preinstall compilation logged default -O2. Fresh original source archive at `/private/tmp/sceptre-cl3-reference-source-627e7e4` generated `/private/tmp/sceptre-cl3-oracle-627e7e4` with debug=false. All four thread variables were one and source ondisc was 1.3.5. Stale existing oracles were rejected rather than weakening environment verification. Strict provenance/toolchain/RNG/dependency/threads checks pass in the new oracle-backed contracts.

- Focused: 16 tests, 195 passing expectations, zero failures/errors/warnings/skips.
- Existing contracts/native/assignment/parallel/auxiliary: 36 tests, 811 passing expectations, zero failures/errors/warnings/skips.
- Combined: **52 tests / 1,006 expectations**. Parse checks and git diff whitespace checks pass. Static current-source/test search finds no retired private symbols; current metric selector updates belong to the root.

This is focused qualification of CL3, not final acceptance of the complete refactor. Parent source review, integration, full immutable package/CI/consumer/resource gates and the coordinated native rename remain. No resource timings or performance claims were made in this lane. The granted build/test window was explicitly released immediately after both processes exited.

The initial broad write was automatically rejected before execution. After rereading the exact signed plan, smaller module-scoped requests were independently reviewed and approved with explicit invariant evidence; no alternate mechanism was used to bypass the rejection. Exact scope/reason is recorded in `/private/tmp/sceptre-cl3-approval-receipts.md`; the machine tool transcript is authoritative.
