# Execution and mixture cleanup proposal — source c7d6250

Observed 2026-10-16. **Proposal only; no implementation.** The dedicated worktree
`/Users/ekatsevi/.codex/worktrees/refactor-execution-cleanup/sceptre` is clean on
`codex/refactor-execution-cleanup` at
`c7d62502757f2a17741e2ef18b9e4d006bf063dc`. No edits or stash exist. The companion
`/private/tmp/sceptre-execution-cleanup-inventory-c7d6250.json` records all 19 named
functions (15 top-level and four local workers), actual source callers, current
formal counts and hashes of the five inspected modules. Anonymous error/warning
handlers in the mixture module were also inspected.

This audit reread all five modules, their S4/assignment/context call boundaries,
the 11-script consumer inventory and relevant test sources. Source parsing was
executed; no package workflow, pipeline, benchmark or test was executed for this
read-only audit. Existing qualification of source 28df390 remains historical
and does not qualify any proposed cleanup.

## Actual call graph and boundaries

- `R/s4_analysis_functs_2.R:175,278,381` call the private S4 analysis adapter in
  `R/analysis.R:3`. It is the only caller of the two drivers, through string
  `do.call` at analysis.R:61,63. The adapter creates the permutation bank,
  performs initial GC and assembles the complete result family.
- Permutation driver `R/association_execution.R:43` traverses unique responses,
  with local block worker at :98. It calls the shared-fit and NT group-fit
  workers at :146 and :160. One response is read per outer iteration; shared
  response fits and score pieces are reused across requested guide groups.
- CRT driver at :223 has a response-fit prepass worker at :277 and guide-block
  worker at :329. Calls at :356 and :376 each create an assignment fit and full
  B1+B2+B3 bank **once per guide**, before looping over its responses. NT-only
  discovery requires a different treatment+NT fit for each guide.
- `R/assign_grna_functions.R:27` calls mixture orchestration at
  `R/mixture_model_functs.R:1`; its block worker at :54 calls the numeric mixture
  function at :79. The latter is only called with `use_glm=TRUE` (:84).
- `get_id_from_idx` is used by three association loops (:108,286,338) and mixture
  (:69). Its apparent ID lookup conceals progress/file writes and GC.
- All 11 pinned pipeline script rows in `tools/refactor/consumer-inventory.json`
  were checked: none call the five modules' numeric drivers/workers/helpers.
  Pipeline private compatibility obligations lie in S4/QC/calibration/result
  processors, outside these proposed private renames. The book contract uses
  exported workflows. This does not establish that arbitrary third-party code
  never calls undocumented internals.

## Decision register

Every disposition below requires inclusion in the reviewed cleanup addendum
before implementation. No compatibility alias should be added solely for a
private name with no actual inspected consumer.

| ID | Evidence and problem | Proposed disposition and boundary | Behavior/risk gate |
| --- | --- | --- | --- |
| E1 | analysis.R:31–63 passes a context plus the same store/design/pairs/cache/index construction inputs; string dispatch hides actual calls. Permutation has 24 formals, CRT 23. | Call `run_permutation_analysis(context, execution, synthetic_idxs)` or `run_crt_analysis(context, execution)` directly. `execution` is one ordinary list of output/progress/parallel/log settings, not a class or framework. Context owns data/design/index/pairs/test settings/cache. | Identical scientific outputs, IDs/order/factors/decisions, messages, warnings, RNG, cache updates and resource scope. No raw assignment/cell arguments duplicated into drivers. |
| E2 | Both drivers accept `context=NULL` and reconstruct an index (:78–82,:257–261), although the only caller always supplies it. | Require the prepared context; delete fallback construction and now-unused arguments. Both `n_nonzero_*_thresh` formals are entirely unread. | Migration atomically updates only known caller, drivers and tests/tool selectors. Do not move pair eligibility back into inference. |
| E3 | NT permutation worker :58–62 and NT CRT worker :203–209 reconstruct partial indices from historical raw arguments. Both actual calls already supply the real index. CRT :144–147/:225 additionally rematches rows or composes cells when explicit values were omitted. | Delete nullable index/response-row/source-cell fallbacks. Driver resolves rows once, passes explicit row indices and fit view/source cells. Remove sole call to `legacy_fit_source_cells`; coordinator can then remove its definition in analysis_context.R:165. | Do not change ordering/duplicates in views, source-to-fit mapping, observed one-based versus synthetic zero-based native coordinates, or callers' owned index vectors. |
| E4 | Existing names `*_in_memory`, `*_v2`, `*_factored_out` describe history or an implementation accident. Backed data already use these paths. | Rename private drivers/workers by traversal and fit responsibility; see proposed interfaces below. Rename S4 boundary `analyze_association_pairs`; drop its redundant `calibration_check` formal and use context's analysis role. | Coordinator must update the three S4 caller names/arguments. Exported S4 functions and formals stay identical. No alias layer. |
| E5 | get_id_from_idx :2–38 combines lookup, display/file output and 200-step GC; mixture imports this association-named helper. | Put plain `as.character(ids[position])` lookup in each caller. Add explicitly side-effecting `report_worker_progress(...)` and named periodic worker maintenance in execution.R, shared by all four loops. | Preserve order: progress before periodic GC before row load/fit; print at 1 and multiples of 5; GC at multiples of 200 even if progress disabled. Keep exact bytes/newlines of message versus write, labels, filenames and file create/append lifecycle. |
| E6 | Both drivers call lookup_response_fit twice on cache hits (:129–130,:300–301). | Store one lookup in a local variable, branch on that value. | Existing read is pure; keep response read before lookup in CRT prepass, not an opportunistic resource change. Preserve cache provenance/insertion order and fit view. |
| E7 | Current design subsetting and NT-only single-column dropping are duplicated in the two drivers (:90–95,:269–274). | Adopt the bounded matrix-preserving `drop=FALSE` correction before or atomically with CL3. Coordinator makes common fit-view/design preparation one context-boundary responsibility, retaining one design in its actual common-view order and one named view. Group-fit slicing also preserves matrix dimensions. The saved full-source design and cache invalidation stay unchanged. | Existing source/audit evidence records the old intercept-only failure; it is not the reviewed target. Add independent intercept-only NT calibration/discovery tests for both traversal modes, with expected views/calculations independent of the new context helper. Keep frozen ordinary-case fixtures unchanged. No all-guide cached matrices/views. |
| E8 | obtain_em_assignments:143 has optional `use_glm`; only call uses TRUE. With FALSE, g_mus_pert0 is undefined before native EM (:166). | Separate decision: remove unused private false mode and formal, leaving the existing Poisson initialization unconditional. Do not implement a new no-GLM mode. | First verify complete source/test callers (done for this snapshot). Run actual mixture workflows and RNG assertions. Private direct unsupported false-mode behavior is not a public feature; record deletion plainly. |
| E9 | Mixture :179–180 falls back to raw-count threshold for low counts/nonconverged EM. This can superficially resemble compatibility fallback. | **Retain**; it is an active scientific behavior, as are expected Poisson-warning muffling, formula diagnostics and fixed starting guesses. | Exercise EM-success, low-count fallback and nonconvergence paths; identical assignment sets/order and conditions. Avoid numerical algorithm changes. |
| E10 | execution.R centralizes real serial/fork dispatch, block binding and progress-file lifecycle. One-line binder is small but used four times and owns schema/factor behavior. | Keep these focused responsibilities; no blanket wrapper deletion. Mixture may reuse file/dispatch helpers only if exact announcement text, flattening and return shape stay explicit. | `rbindlist(fill=TRUE)`, `mclapply` block order/core count and two-phase CRT log append must remain. No new scheduler abstraction or generic callback registry. |
| E11 | No named internal function in the five modules currently has @noRd contracts; stale skipped numbered comments and “combine and sort” at :206 misdescribe code (there is no sort). | Document each retained named internal function, including four local workers. Explain input coordinate/view, returned values, active cache/bank scope and side effects; replace stale numbered/history comments with purpose comments. | Contracts must match code; do not describe workers as in-memory-only, pure if they read storage, or normalize/sort rows to make a test pass. |

## Proposed interfaces for review

Names and field ownership below are concrete proposals; they have not been
implemented or agreed by a senior reviewer. Use the current small ordinary
`context$test_options` (B1/B2/B3, side, approximation), not a giant catch-all
options object. `output_amount` belongs to ordinary execution options.

```r
analyze_association_pairs(
    sceptre_object, response_grna_group_pairs, analysis_type,
    output_amount, print_progress, parallel, n_processors, log_dir)

run_permutation_analysis(context, execution, synthetic_idxs)
run_crt_analysis(context, execution)

# Called with already aligned response and score arrays; no S4 facade.
test_permutation_shared_fit(
    expression, score, treatment_lookup, groups,
    synthetic_idxs, test_options, output_amount)

# Prepared index is the one owner of coordinate mapping. Each active group
# materializes its treatment+NT view only when needed, not for all pairs upfront.
test_permutation_group_fit(
    expression, design, index, groups,
    synthetic_idxs, test_options, output_amount)

# Explicit row indices/source cells; response IDs retained solely for cache keys.
# Creates one assignment fit/full bank before traversing the requested responses.
test_crt_shared_fit(
    response_store, response_rows, response_ids, design, source_cells,
    treatment, response_fits, test_options, output_amount)

test_crt_group_fit(
    response_store, response_rows, design, fit_view,
    test_options, output_amount)

# ID lookup happens visibly in the caller. This helper only reports progress.
report_worker_progress(
    position, id, total, print_progress, parallel, log_file,
    feature = "response", prefix = "Analyzing pairs containing ")
collect_worker_garbage(position) # fixed 200-item cadence, no user-facing option
```

Context preparation remains one S4-to-ordinary-data boundary. Proposed context
fields are one response store, index, one common-view design, one common view,
eligible pairs, role, small test options and existing response fits. Preserve
only fields that actual code needs. A treatment lookup closure is prepared once
from that index; coordinator-proposed rename
`prepared_index_factory` → `make_treatment_index_lookup` must update both drivers
and definition together. This is not a registry or extension mechanism.

Workers remain four small scientific cases rather than one configurable engine:
permutation/shared fit, permutation/group fit, CRT/shared fit, CRT/group fit.
Storage/fitting/scoring/native inference/cache ownership stays in the existing
modules. The native sampler/inference APIs are untouched by this proposal.

## Migration sequence and stop points

1. Obtain senior review of E1–E11 and context/design ownership. Confirm deliberate
   retained behaviors, the adopted bounded intercept-only matrix correction,
   dead-mode deletion and protected public/consumer boundary. The original E7
   proposal to preserve a vector/error is superseded by that correction.
2. Single coordinated source unit: rename S4 adapter and its three callers;
   require context in both renamed drivers; direct calls; remove duplicate fields
   and threshold arguments; one common-view design owner. Update private tool
   selectors and current source-tracing docs, retaining historical reports.
3. Second unit: four explicit worker interfaces, remove index/row/cell fallbacks,
   remove unused NT raw assignment inputs, coordinate factory rename and deletion
   of sole legacy_fit_source_cells. Keep traversal/bank generation scope fixed.
4. Third unit: progress/maintenance separation and mixture dependency cleanup;
   separately approved removal of unused use_glm false mode. Keep real fallback.
5. Contracts/names/comments reviewed alongside each unit. Refresh call graph,
   deleted-path register and static metrics after stable interfaces. A reduced
   formal count alone is not proof of maintainability.

Each unit must be internally coherent and independently revertible. No aliases
or temporary fallback branches should remain at completion. If context work and
worker signatures cannot land coherently separately, combine those changes in
one reviewed unit rather than leave an invented compatibility bridge.

## Tests and resource/RNG gates

Existing test sources inspected (not run here):

- test-refactor-contracts.R:224 exact exported API/S4 facade; :229 supported
  workflow scientific values with 1e-8 numeric tolerance and exact identity/
  decisions; :258 ordered index mapping; :267 cache behavior; :282 old saved
  objects; :289 degenerate behavior; :298 serial/two-worker behavior; :316 backed
  permutations and independent attachment.
- test-refactor-indexing.R:19,56,84,96,140 covers noncontiguous cells, native bases,
  NT overlap/order, composed contexts and reordered requests.
- test-refactor-cache.R:1,22,50 covers view-valid cache reuse/invalidation.
- test-refactor-native.R:52,83,101 covers sampler laws and exact pinned draws.
- test-assign_grnas.R:424 active mixture assignment; :493 enriched formula error;
  test-parallel_mac.R:57 mixture fork execution; test-aux_functions.R:160 seeded
  starting guesses.

Targeted additions should test behavior missing from these sources, not mirror
new implementation lines:

1. Capture serial progress and raw parallel progress-file bytes at positions
   1/5 and around 200 for association and mixture; progress disabled creates no
   progress file and still performs periodic worker maintenance. The cadence
   counter restarts per worker/block, and CRT response prepass versus guide phase
   are separate. Instrument GC only in this narrowly scoped contract test.
2. Trace fitter/sampler calls in actual small public workflows to verify one
   response fit for a shared view, per-guide NT response fits, one CRT assignment
   fit and full-B bank per guide. Compare true output against the pinned oracle,
   not mocks alone; trace call order/draw budgets rather than implementation names.
3. Exercise mixture success plus active count/nonconvergence fallback and verify
   caller RNG preservation and assignment ordering, serial and two-worker cases.
4. Rerun existing scientific/index/cache/native/saved-object contracts after
   each coherent unit. Exact IDs, request/pod/order, schema/factors, stage/QC/
   significance decisions, warnings/errors/messages and RNG must be compared;
   only numeric scientific fields use 1e-8 tolerance. Do not sort payload rows or
   factor levels to conceal changes. Existing cache keyed-comparison behavior is
   explicit and does not authorize sorting consumer outputs.
5. Rerun touched real consumer paths after installing an immutable source build:
   the 30 actual-script paired stages and affected uncached book cells; full
   Nextflow only if affected assembly/orchestration behavior requires it. Earlier
   28df390 evidence is historical. Preserve the baseline completion-order caveat
   and known CellRanger barcode failure rather than mask them.
6. Resource qualification uses matched -O2 reference/candidate, fixed threads,
   source ondisc 1.3.5 and unchanged build dependencies. Watch peak memory and
   elapsed time on the existing sparse/backed fixtures: no response matrix
   materialization, no all-pair/all-guide view bank, no response refits inside
   shared permutation group loops, no CRT bank creation inside response loops,
   and no lost GC/row-read behavior. Do not claim performance gains from counts.

The highest-risk changes are common-view design ownership (including the adopted
NT intercept-only matrix correction and supported overlap/order contracts), worker
treatment/source coordinates, bank
scope/draw order, progress-file lifecycle and cache insertion/merging order.
Plain renames and unread-argument deletion are lower risk but still require the
same public API and source-call readback. No new extension is part of this proposal. The intercept-only correction is a
named bounded behavior change requiring its own independent regression; historical
error evidence must remain distinguished from the reviewed target.
