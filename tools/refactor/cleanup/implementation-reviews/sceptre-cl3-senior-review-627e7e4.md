# Independent CL3 source review

Reviewer: rf0_senior_review. Date: 2026-10-16.

**Outcome: no remaining P0, P1 or material P2 runtime/source finding in this bounded CL3 implementation.** This is not final installed, consumer, CI or resource acceptance. The reported 195 focused expectations and currently running original-oracle contracts do not substitute for those gates.

Reviewed `627e7e4c619b06359980b9e04fb246faa67c398d` versus `f8f55b1a546dd44261ed0eaa205c622beb09dcea` in `/Users/ekatsevi/.codex/worktrees/refactor-rf2-storage/sceptre`; HEAD and the clean working tree were checked. Every added/modified source/test below matched its immutable Git blob. The deleted adapter file was read from the parent. Frozen fixture paths have no diff. No production edits, test reruns or benchmarks were performed by this reviewer.

## E1–E11 and coordinate ownership

The three S4 callers now use a direct `analyze_association_pairs` boundary and pass no duplicated calibration flag. That boundary constructs one context and one small execution list, creates the same shared permutation bank, performs the initial GC and delegates to one of two explicit traversals. Drivers require prepared data; nullable context/raw argument reconstruction and unread thresholds are removed rather than relocated. No public formals or S4 slots are changed.

`prepare_analysis_context` is the sole owner of the common fit-cell design slice, retaining `drop=FALSE`. The object keeps its full source-row design and compact saved fit cache. Common source cells are all analysis cells for complement/discovery, or the ordered NT pool for NT calibration. Each NT discovery group alone composes treatment-first plus NT positions and selects its own matrix design/response view. No all-guide bank of designs/views or full assay is created. Tests supply independent source/treatment expectations and an intercept score, cover calibration/discovery under both mechanisms and verify object serialization. The bounded intercept-only correction is correctly separated from historical error evidence; the immutable edge fixture is retained and only that named expectation is superseded.

`encode_nt_fit_assignments` preserves supported disjoint low-MOI saved NT partitions, pool order and names with explicit membership-to-pool matching. Actual QC callers are updated; the one-call forwarding alias is deleted. Unsupported high-MOI NT controls remain prohibited by public validation. Overlap/union behavior is tested independently rather than retaining cumulative out-of-pool positions as compatibility. Native observed positions remain one-based at the actual Rcpp boundary; synthetic bank coordinates remain zero-based. Deleting the integer passthrough does not change valid coordinate types.

## Traversal, cache, bank and effects

`run_permutation_analysis` retains first-request response order and within-response pair order. `run_crt_analysis` retains its shared-response prepass followed by first-request guide traversal and within-guide response order. Partitioning source-row/group integers rather than label vectors is equivalent for validated unique labels: seeded shuffling operates on the same ordered positions, and the partitioner only calls sample when length is at least twice a positive worker count, avoiding numeric singleton sampling semantics. Factor labels/levels remain in the original pair tables used to construct results. New tests compare label-versus-integer block positions and caller RNG.

Shared response fits are read once from the cache per lookup point, computed only on misses and merged in existing block/traversal order. Both permutation outer loops and CRT prepasses still read responses before cache lookup, including warm hits. NT discovery response fits stay per pair and uncached. The cache payload schema, provenance/invalidation and source-row design contract are unchanged; its prior misleading design-order comment is corrected.

The four workers remain distinct scientific cases. Shared permutation scores are computed once per response; guide-specific NT permutation fits remain within that response's group loop. Each CRT guide fits assignment probabilities and draws one complete checked B1+B2+B3 bank before response traversal; group-specific responses use the correct current fit view. Native stage budgets, approximation, output-distribution flag, side and use_all_cells are passed unchanged. No new lazy sampler, partial bank or fitted-vector cache is introduced. Actual-workflow instrumentation tests enforce fit/assignment/read counts, full bank budgets, warm behavior and RNG rather than only mocking return values.

ID lookup is now visibly pure. All three association loops and mixture call reporting, then periodic GC, then storage/fitting in that order. Cadences remain first/every-five progress and every-200 maintenance independent of reporting. Existing log helpers retain create/truncate versus CRT append behavior, file names, announcements, extra write newline and serial/fork block order. Helper byte/cadence contracts and a real 402-guide mixture block test cover the relevant mechanics. The latter tests deterministic partitioned names, assignment equality, log bytes and RNG. Source inspection also confirms unchanged CRT two-phase log lifecycle; final serial/fork workflow qualification remains necessary.

Mixture initialization remains the active Poisson path. Removing the unused use_glm false branch does not replace it with a new mode. Native EM, expected Poisson-warning muffling, contextual formula error, fixed starting guesses and raw-count fallback on low counts or unsuccessful EM are retained. Tests execute actual EM success and actual nonconvergence, plus the low-count fallback without an EM call. Fitter/score arithmetic is unchanged; unused coefficient/reduced-score modes and the three precomputation forwarders are removed with all actual current R/test callers migrated.

## Source clarity and integration dependency

The change simplifies ownership: one ordinary context, explicit execution effects, direct driver dispatch and four purpose-specific workers, with contracts describing real coordinates, scope, storage/cache effects and retained fallback. It removes historical/private aliases and reconstruction paths rather than introducing another forwarding layer. Search found no deleted-name executable use in current R/src/tests. No gratuitous new fallback, broad framework or scientific policy was found.

The worker branch still contains old selectors in `tools/refactor/measure.py:98-100`. This is **resolved on integration** by independently inspected commit `1eddbe8aa223e5c5cd7d70fbbe062c60bf89574f`: it includes the new adapter, both drivers, four workers and eventual native name while retaining old names only for historical source measurements. Integrate that dependency before regenerating current metrics; do not present worker-branch trouble_spots output as complete. No production alias is needed. Historical measured JSON/spec discussions remain historical and should not be rewritten as current evidence.

Final acceptance still needs the coordinator's exact immutable original/frozen/saved/cache/backed/serial-fork tests and affected consumers, plus matched release resource measurements of the changed common-view/worker scope. Source inspection establishes no added full assay/all-guide materialization or obvious lost reuse; it does not measure RSS/time or prove pristine source.

## Reviewed plan SHA256

- `tools/refactor/CLEANUP_PLAN.md`: `fae06927ca2ec4ae0c49a50a858f61b791ee323320471d049eeede9c6160c7b9`
- `tools/refactor/cleanup/execution-proposal.md`: `95942e0fe2a4eec741cd2cac2c63ade0edda800228ce0a20f62762a67ed9ba4e`
- `tools/refactor/cleanup/symbol-dispositions.md`: `7d382d350862ccb47112f77be091e4580421cbb16b1e5715d529d417f05a4758`
- `tools/refactor/cleanup/senior-review-amendment.md`: `614bbc3e617421cb4df98ae47fd1b1ed49f3db3ab2493adefea7672d2fa36476`

## Changed-file SHA256

Deleted-file hash is the parent's content.

| File | Change | SHA256 |
|---|---|---|
| `R/analysis.R` | M | `05de77ba0b217b64d755a10f8eceee60bc472f029dbb1b0ea4e1e0a0a93b2497` |
| `R/analysis_context.R` | M | `c24034bbc966a1bdd7f49df99dee4f58973279a0e324e7baeee2a839f1a456b7` |
| `R/association_execution.R` | M | `8915709e29869860bae5b6562b965318063a33c1dfdf161369d62e07ca022451` |
| `R/association_workers.R` | M | `dc9f8e5722463734e717c49dfa91af8035deb93aaeb179e2bb970c336dcdc58b` |
| `R/cellwise_qc_functs.R` | M | `a650fa3cbc82610d774c730f15565c7ce0a774316b47791818db2452a0b422db` |
| `R/execution.R` | M | `679d96e3f20e8be9178cdc4de1e3de3dd7503d941a49db8687c93fe4c4981235` |
| `R/fit_cache.R` | M | `9decdf97b39ec9da539bf806cd4064ff033c60e0109254064cf0cb6173ce22a1` |
| `R/fitting.R` | M | `0b6935e9439ff9e32f4699ad1968c0fdd8180d7c87f955b255b226fe77208993` |
| `R/mixture_model_functs.R` | M | `d3d05730b8912169063dcfb0ffea94056ab188fd276b355470537b5ec5605ff8` |
| `R/precomputation_functions.R` | D | `ab5923e482adf963da55b9a1bd23459239763c940d5891500aaad1f1088ee5ec` |
| `R/s4_analysis_functs_2.R` | M | `1e11d0b35e75928d12b7790c65df01c46d317afd24a9d771d2d213c39b9b626a` |
| `R/score_preparation.R` | M | `9d7d9088c452bb7a28e71290184916aec990734464fe3d905c3d8a1c08b64bb1` |
| `tests/testthat/test-refactor-context-execution.R` | A | `36d2630b1b610751b0f40e3ffa87f6fb9805ca96541bdd702847094ac1defb48` |
| `tests/testthat/test-refactor-contracts.R` | M | `0e9e18bf089020ac0edef76aaf2a59f85d1baf0f7468d86f5c0e300e8850d36b` |
| `tests/testthat/test-refactor-fitting.R` | M | `c3fd9a0e679cd1dff1b2af3d6efd676d087dfff5c95bd4a4b676a42d342db056` |
| `tests/testthat/test-refactor-indexing.R` | M | `0f4eee47712b0383e26759f8a9fc005306e22e08a5f794ac9e483037a991c64e` |
