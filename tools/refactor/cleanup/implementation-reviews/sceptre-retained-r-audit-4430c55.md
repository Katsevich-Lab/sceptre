# Retained R contract and cleanup closure audit

Target: integration `4430c55a3f05780f7ac6ec3bd380b1ecc198abe2`. Initial readback HEAD was `b1a4c78255db7e97ed3e090b17d4aa941a21f0f5`, whose R source was identical to target. Root subsequently committed the six reviewed documentation-only R amendments and numerical-owner comments as `750b1b5da81de3fb3189f2ca95368c4beaa9c618`. Every retained R hash below was checked against that exact commit. Native comment review is bounded separately below. No production edit, R process, build, test, benchmark or external action was performed by this reviewer.

**Outcome: no remaining P0/P1/material P2 in the retained R source/contract scope after the reviewed local documentation amendments. This is not final source or resource acceptance.** Native CL4, the explicitly active B6 aggregation work, the distinct dense-entry repair and their integration qualification remain open.

## Findings and actual closure

1. At the immutable target, `R/cellwise_qc_functs.R:1` and `R/s4_analysis_functs_2.R:89,160,620` lacked useful contracts despite owning cell selection, two actual pipeline calibration phases and automatic workflow advancement. The local comments now describe source-cell order, overlapping reason counts, QC/cache responsibilities, prerequisite/default side effects, validation/sampling order, raw versus finalized pipeline results and successful facade/cache writeback. They agree with actual bodies and callers. The quantile closure shares its enclosing contract; adding repetitive standalone roxygen is unnecessary.
2. `R/getters_and_setters.R:76` promised a sparse matrix while its intentionally preserved QC subsetting can drop to a logical vector. The local return-prose/Rd amendment discloses that existing behavior; it does not silently fix the preserved dimension drop.
3. The permutation-bank chooser contract previously deferred its essential meaning to an unspecified view. The local `R/preprocessing_functions.R` detail now explains its four real view/size branches and the hybrid guarantee **only for m <= k <= M**, full-bank creation/reuse and treatment-prefix consumption. This is the R owner's actual policy boundary, not a new generalized sampler.
4. The local Parse loader comment now states cell-first file coordinates, owned zero-based conversion, feature-by-cell output and delegation to the complete CSC/duplicate-rejection constructor. `process_calibration_result` now correctly calls its skipped Bonferroni collapse non-calibration (power/discovery), rather than discovery-only.
5. The architecture guide originally named a retired integration function/obsolete delegates and overstated mutable native-counter inputs. The current local corrections agree with retained R ownership. A final whole-guide audit is still required after native names and B6/dense behavior settle; this interim pass does not certify that future final text.

These are source/prose findings, not conclusions from expectation counts. The historical `legacy-contract-checkpoint.json` is a pinned checkpoint; its old noRd counts/names do not substitute for final source coverage.

## Scope and coverage

The original inventory contains **34 handwritten R files**. Current source has **32**: `precomputation_functions.R` and `test_functs.R` were retired. `R/RcppExports.R` is generated, not a handwritten module. Dataset/package/S4-schema files and three S4 display/plot methods were included; local named closures were inspected under explicit enclosing or separate contracts. Anonymous callbacks/condition handlers are not falsely counted as missing named owners.

Current contract text, actual R references and the pinned consumer inventory were checked against the earlier independently reviewed source units. The attached manifest records each retained file's exact target/working hash and byte matches to prior reviewed unit checkpoints. This allowed reuse of actual earlier body reviews while checking the final names/comments/callers. Hash matching and regex navigation are not semantic test coverage or proof of absence of bugs.

| Files | Responsibilities checked |
| --- | --- |
| analysis.R; analysis_context.R | One required context; source/analysis/NT/fit views, order and repetitions; request identity; lean treatment lookup; no per-pair full-map or all-guide view bank. |
| association_execution.R; association_workers.R; execution.R | Distinct response versus guide traversal; cache-hit reads; one guide-bank scope; cache miss/insert order; progress then GC then read; serial/fork order/log lifetime. |
| fitting.R; score_preparation.R; fit_cache.R | Compact nuisance records, matrix orientation/eigen order; actual fit-view validity; full-design/control/assignment/cell invalidation and successful commit. |
| assign_grna_functions.R; cellwise_qc_functs.R; mixture_model_functs.R; pairwise_qc_functs.R; pair_qc.R | Real maximum explicit interface, assignment pooling/QC views; union/duplicate/stored-position semantics; shared eligibility; starting-guess order/restored RNG; active mixture fallback. |
| s4_analysis_functs_1.R; s4_analysis_functs_2.R; workflow_state.R | Public facade/formals; real pipeline split stages; validation/default and workflow transition order; successful slot/cache commit; sibling final roles. |
| import_functs.R; storage.R; ondisc_functs.R; getters_and_setters.R | Owned CSC completeness/ranges/repeated-coordinate rejection; owned-versus-external parser boundary; source order, native/backed coordinates, static copy ownership; connection cleanup and reopen identity; preserved getter edge. |
| neg_control_functions.R; preprocessing_functions.R; aux_functions.R | Post-QC authority/default rejection; scalar budgets/bank bounds; calibration sampling/group labels; actual design and metadata joins; deterministic partition and log contracts. |
| result_processing.R | Native amount/schema conversion; failed rows; real guide integration/Bonferroni roles; explicit alpha outside amount-2/3 mask; complete-family pipeline adapters. |
| plotting_functions.R; qq_plot_helpers.R; pair_constructor_functs.R | Forwarded appearance flags; stored significance; per-panel QQ grouping/RNG; preserved display conventions; wide/inclusive cis rule and owned trans exclusion. |
| s4_helpers.R; s4_class.R; data.R; sceptre.R | Real generic methods/lookup boundary, persisted schema and data/package documentation; no fabricated runtime contracts for data-only files. |

## Explicit remaining handoffs, not forgotten R cleanup

- **CL4:** migrate four `run_low_level_test_full_v4` worker calls and two `crt_index_sampler_fast` calls in `association_workers.R`; rename the calibration constructor/caller and its native sampling/increment calls in `neg_control_functions.R`/`s4_analysis_functs_2.R`; regenerate both binding files, remove retired registrations and migrate executable harness/test lookups atomically. Current old spellings are expected before that unit, not final closure. `increment_vector` has no remaining handwritten R caller; its native retirement is still pending.
- **B6:** `collapse_grna_matrix` still has its full-dense implementation/temporary rewrite comment at this target. Existing summation order, metadata repetition and numerical state remain the characterization constraints; sparse aggregation and meaningful resource qualification are actively assigned. Do not claim the final R source has no avoidable dense realization yet.
- **Distinct dense repair:** `bc53bc1` was separately source-reviewed and fixes the existing symmetric/near-symmetric entry loss with supported general-before-T coercion. It is not part of this immutable target, and its installed qualification/dense-copy measurement remain pending.
- **Final acceptance:** reread integrated R/native contracts and all handwritten implementations/headers, whole ARCHITECTURE.md, exact API/56 slots/old-RDS/RNG/order/cache gates and final resources on one immutable combined source. Unchanged scientific histogram labels, Bonferroni guide/target interpretation, plot sampling caps, backed CRT and PR233 remain declared separate decisions; none is silently introduced as cleanup.

No additional material wrapper/fallback/coupling cleanup is justified by this pass. Retained storage/QC responsibility boundaries, the lean lookup closure and actual private pipeline adapters have real consumers. Their small size alone is not a reason to delete them.

## Evidence binding

Manifest: `/private/tmp/sceptre-retained-r-source-manifest-4430c55.json`, SHA256 `aa1e93a90b5a99ab2b4ace4981e549a1aa18410d552f8c666956538fe41a8aef`.

- `R/analysis.R`: `05de77ba0b217b64d755a10f8eceee60bc472f029dbb1b0ea4e1e0a0a93b2497`
- `R/analysis_context.R`: `c24034bbc966a1bdd7f49df99dee4f58973279a0e324e7baeee2a839f1a456b7`
- `R/assign_grna_functions.R`: `d9ce47de1d9ba371ef366deda9c26744d6616d220273b3e19c0f32c688a5f48c`
- `R/association_execution.R`: `8915709e29869860bae5b6562b965318063a33c1dfdf161369d62e07ca022451`
- `R/association_workers.R`: `dc9f8e5722463734e717c49dfa91af8035deb93aaeb179e2bb970c336dcdc58b`
- `R/aux_functions.R`: `e19d9dbc78f091f0e71347aae299c9a0228e8f6526210a01c5bb744abc653dc6`
- `R/cellwise_qc_functs.R`: `97f7dd7c715f54472c5eccfe1e7b7d3883aea68bb6c5cd2c4cf6fbbb827862ff` (documentation amendment committed in 750b1b5)
- `R/check_functions.R`: `3b021705dcf19128b263e871fe65ee75bd45ed0a8630b238b859367560731f34`
- `R/data.R`: `8a35c3f2206b3424701158e5bb6315d125d10cc5535820406c6c4d69432cc4a6`
- `R/execution.R`: `679d96e3f20e8be9178cdc4de1e3de3dd7503d941a49db8687c93fe4c4981235`
- `R/fit_cache.R`: `9decdf97b39ec9da539bf806cd4064ff033c60e0109254064cf0cb6173ce22a1`
- `R/fitting.R`: `0b6935e9439ff9e32f4699ad1968c0fdd8180d7c87f955b255b226fe77208993`
- `R/getters_and_setters.R`: `26647c51ce4c9b855c08e9d3bd01f5308239a384555ef26a3134259c1aa59f41` (documentation amendment committed in 750b1b5)
- `R/import_functs.R`: `b2fa50e9c82205457ac13d65fd7b528fef47739775706e80416d22db609054c7` (documentation amendment committed in 750b1b5)
- `R/mixture_model_functs.R`: `d3d05730b8912169063dcfb0ffea94056ab188fd276b355470537b5ec5605ff8`
- `R/neg_control_functions.R`: `f4fe402cc769cb8d823cf31f18a1b86c2c2a57f664a78aa20c19fce8cb6717dd`
- `R/ondisc_functs.R`: `5c6579aac2477259e2209cc0da68ccce839971abab96df08015f83decc19704d`
- `R/pair_constructor_functs.R`: `95135d4871d73b90e53bba1f1f06cb3ae4210ddadc17155c2198ef7c685b7971`
- `R/pair_qc.R`: `dc77c0f53390416e4ba3c881167fc73485bc578f64871ac5b7a7142467629be0`
- `R/pairwise_qc_functs.R`: `27f9ecb236ef0b56dbd2ebc880bc847201bd31b3d15fa1d91ad6eb516cc4a5e8`
- `R/plotting_functions.R`: `66b1547678d785da9d14ed2fd784ac8d5ed42bb88059a73742bb65c8de98bcfa`
- `R/preprocessing_functions.R`: `0658ba6fa9737ade671633bef93cda547e6b693dec6ed69845b6d0ae0806eca1` (documentation amendment committed in 750b1b5)
- `R/qq_plot_helpers.R`: `348fc17cde0864bea749054fa755324d7f1b24dd746395cd6eb8c6b19ec7d568`
- `R/result_processing.R`: `5d94e22022669ba12af53b7dd221d05c3149ca367b857625f1807e8072b21185` (documentation amendment committed in 750b1b5)
- `R/s4_analysis_functs_1.R`: `52df97f69fe6f5a9b9876ac4bbccb92fab7c03df88b1ad31d8d8ce1c73b8d3d1`
- `R/s4_analysis_functs_2.R`: `026e88e022d7a93e09ddb34d813524d121cdeafab3851d2bc86f45f120c55868` (documentation amendment committed in 750b1b5)
- `R/s4_class.R`: `3f66cd9c563e6e51ef1b0e838831b404bc52f08545494f48b11ea026f6539998`
- `R/s4_helpers.R`: `29a31262d8cae42b7dc48fd29d42897c58ba23a9a098be1c8db40a07a5d5d770`
- `R/sceptre.R`: `a38a62c6bef1af6a1b52ea60edabcde900eb71de6e23b1cbca4f25a941a62ecf`
- `R/score_preparation.R`: `9d7d9088c452bb7a28e71290184916aec990734464fe3d905c3d8a1c08b64bb1`
- `R/storage.R`: `3755dea9c81629b0140b03230a09d6359e57169ea02580ce4f19c1da9b4df06e`
- `R/workflow_state.R`: `8f3d2e01033ea89e7f46a8faec76f3f902fe374b956e5955050584cd897443ca`

## Numerical-owner comment readback at 750b1b5

Read both complete native files and the comment-only diff, with no numerical execution. The fixed-Poisson-mean NB variance parameterization is correct. The information expression is minus the derivative of the score, so the stated plus Newton step is accurate. The MM expression matches the Pearson residual sum to supplied residual degrees of freedom with the implemented denominator. MM starts at the original pilot; its warning restores that pilot, while an exception before its result commit preserves the preceding MLE estimate and method 1. These comments describe the actual bounded typed catches and partial state, not a new fallback.

The mixture comments correctly describe baseline/additive Poisson means, the density floor solely for likelihood, algebraic posterior log odds, posterior/weight flip above 1/2, deliberately preserved log-ratio update used additively, converged-start competition and strict first-winner tie rule. The R caller's low-count/nonconvergence thresholding remains separate. One minor stale inline phrase was sent to root: “not NA or uniformly 0” should say any nonfinite posterior or all posteriors <=1e-100. That is the exact degeneracy guard and differs from the likelihood density floor. It requires a comment edit only; it is not a material runtime or design finding.

This readback covers these two owning mathematical comments, not all native files/headers. Root's noncomment-token identity artifact was inspected; its new hashes match this commit:

- `src/compute_nb_size_param_functions.cpp`: `69eb74805365148f5a57c8ae1a1801732c4e6a6c6d3ae71ee0d1c21ea1c70628`
- `src/mixture_functs.cpp`: `216af4e9cb91f3aaba547df81661b61aad4351c2281c191607eea87720c635dc`

## Mixture wording closure

Exact commit `88a613afe3dce3670437b72bd547106769f3660d` resolves the remaining inline comment nit: it now states any nonfinite posterior or none above 1e-100, explicitly separate from the likelihood density floor. The winner comment also correctly limits initialized-record fallback to **no winner**, rather than a later failed start after an earlier winner. Independently checked the exact diff and unchanged noncomment expressions. `src/mixture_functs.cpp` SHA256 is `bfd1badf8121a4873d5c3236ae5eb8a56b5967ccf6a15c19147d3e7614fac248`. No remaining finding from this bounded mathematical-comment review. No numerical test was needed or run for unchanged expressions.
