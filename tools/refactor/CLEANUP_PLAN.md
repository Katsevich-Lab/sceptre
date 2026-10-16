# Source cleanup addendum: inventory and review checkpoint

Status: draft for independent senior review; **further cleanup implementation is
paused**. Source inspected: c7d62502757f2a17741e2ef18b9e4d006bf063dc,
October 16, 2026. This addendum supplements REFACTOR_SPEC.md; it does not
declare the current source pristine or the refactor complete.

## 1. Scope and safe checkpoints

Eugene requires systematic understanding and a reviewed plan before further
cleanup, with genuine elimination of avoidable debt. Public APIs, the 56-slot
object, supported saved objects, scientific methods/defaults, RNG streams,
ordering and resource contracts remain constraints. A test of an undocumented
private helper is not sufficient reason to retain that helper.

The integration checkout is clean at c7d6250. Completed production qualification
predates this cleanup: installed 90664a7 passed 3,522 checks, 112 native component
checks, and the paired package-check flags with the same two omitted-vignette
warnings. Hosted c635263 passed all five lanes. The production expressions of
c7d6250 are unchanged from those qualified builds; its added core contracts
still require text review.

Preserved work is distinguished from accepted implementation:

- The pending canonical sparse-QC/owned-triplet unit is saved on
  codex/refactor-canonical-counts-checkpoint, tree
  8f1b8f4e1340096898d7301b10db4ea365fdc69f, and in a recoverable Git stash.
  It is **not integrated or qualified as a combined unit**.
- Core contracts are integrated in c7d6250; one source/design-order description
  was subsequently found inaccurate and must be corrected.
- The AST-identical 39-function legacy-contract draft is isolated at a2e66f0.
  The storage pointer-name cleanup is isolated at 13ff940.
- The execution cleanup is a proposal only; its dedicated worktree has no edits.
  Earlier resource/consumer reports remain pinned to 28df390, never relabeled.

## 2. Review coverage and evidence limits

[review-coverage.json](cleanup/review-coverage.json) lists every source file,
SHA256 and parsed function location: **34 handwritten R modules, 14 handwritten
C++ implementations and three headers**. Generated Rcpp bindings are classified
separately. The R inventory contains 226 function literals: 158 assigned
top-level functions, three other named-argument owners, three S4 method owners,
and local closures. Counts describe the pinned parser, not comprehension.

| Lane | Modules inspected and caller scope | Detailed evidence |
| --- | --- | --- |
| Core/context/results | All 52 named functions in storage, fitting, score preparation, caches, QC policy, workflow state, analysis context, result processing, old precomputation aliases; all 28 existing contracts checked | [per-function inventory](cleanup/core-inventory.json), [12 findings/dispositions](cleanup/core-review.md) |
| Execution/mixture | S4 analysis adapter, both traversal plans, four workers, execution helpers, mixture assignment; their facade and assignment callers | [19-function inventory](cleanup/execution-inventory.json), [E1–E11 proposal](cleanup/execution-proposal.md) |
| Other R/public boundaries | All other 20 handwritten R modules: validators, imports, pair constructors, assignments/QC, getters, plots/QQ, class/package/data declarations and test-only utilities; generated bindings and three test helpers also inspected | [per-function checklist](cleanup/legacy-inventory.json), [L01–L20 findings](cleanup/legacy-findings.json), [review narrative](cleanup/legacy-review.md) |
| Native | All 14 implementations/three headers, their R callers, numerical/storage/sampling ownership and actual pipeline/book boundaries | Native register below; per-file coverage/hashes in review-coverage.json |

The actual consumer boundary was checked against the 11 pinned pipeline scripts
and selected book source, not all possible external private-namespace users.
Private calibration/QC parts, result processors, pod/file helpers and exported
workflows have identified consumers. The internal numerical drivers, historical
fit aliases and dead native diagnostics do not. Static calls, source reasoning
and earlier runtime evidence are distinguished; unsafe native counterexamples
were **not executed**. This inventory cannot prove absence of every bug.

## 3. Material findings and proposed dispositions

The registers are intentionally broader than the implementation queue. “Fix”
below is a proposal until this plan is reviewed. Normal-case expectations stay
frozen. No methodological repair is hidden in a cleanup.

| ID | Finding, actual behavior and contract | Proposed action and verification |
| --- | --- | --- |
| S1 — P1 | High-MOI calibration in preprocessing_functions.R sums the largest individual NT group sizes for M. Overlap can make M exceed n_cells; Fisher–Yates then indexes outside its source vector. Ten cells/two NT groups each covering 1:8 gives M=16 for a valid union of eight. | First safety unit: bound the existing upper bound by n_cells; check native allocation/draw preconditions before indexing. Preserve M and all draws for ordinary valid cases. Add an overlap-heavy public workflow, preserved safe zero-size samples/banks, and invalid-boundary rejection tests; do not execute the old unsafe kernel. |
| S2 — P2 | Calibration group size is not checked against the actual post-QC nonempty NT assignment groups before combinations/indexing; metadata guide counts can overstate that current set. sample_combinations_v2 narrows a potentially nonfinite estimate when p_hat=0; sample_undercover_pairs_v2 multiplies signed n_groups*n_genes before allocation. check_sn_tail needs at least 100 samples; public B2=4999 is safe. | Apply the per-entrypoint contract table below, including post-QC group loss, with checked wide arithmetic before narrowing. Fail clearly when the current algorithm's candidate-rate/allocation prerequisites are unmet; do not invent a new sampling heuristic. Reject unsupported short native skew-normal inputs rather than read past the vector. Tests must distinguish unsafe inputs from unchanged supported workflows. |
| S3 — P2 | Draw-bank producers allocate a raw heap bank before subsequent vector allocations and only register its external-pointer finalizer at return; an intervening exception can leak the partially built bank. | Use local RAII ownership until successful pointer registration, then transfer ownership once. Keep the exact bank type, contents and existing finalizer lifetime; no pointer class/schema redesign. Review both success and exception paths alongside CL1. |
| I1 — P2 | prepare_index_context decodes persisted NT-local groups to analysis coordinates; storage_pair_counts maps them back to NT-local; the owned sparse kernel redirects through the NT pool per response. This round-trip is not mathematically needed. | Owned sparse QC consumes direct analysis groups; convert R positions once before numeric loops. Keep the same fused all-response summaries, response-sorted scan and request-order restoration. NT totals remain the sum of NT-group counts, including overlap; stored zeros count by position. Retain only the backed ondisc input adapter required by its existing private API. |
| I2 — P2 | Owned group-count/cell-covariate kernels decrement R vectors in place; storage adds defensive copies. The private ondisc kernels can also mutate vectors. | Make sceptre-owned native boundaries use owned C++ coordinate vectors/nonmutating inputs and eliminate the now-unnecessary R copies. Keep copies precisely at ondisc calls that still require them. Test repeated calls and serialized input equality, empty/repeated/ordered groups, stored-zero/fractional semantics and sparse/backed parity. |
| I3 — P2 | The cumulative NT-local encoder is a genuine saved-schema boundary, but a new test enshrines an invalid overlapping NT pool: a length-three pool with positions 3:4. Public high-MOI NT controls are prohibited; supported low-MOI NT assignments are disjoint. | Keep the saved schema and order; replace cumulative construction with explicit membership-to-pool positions and remove its one-call forwarding alias. Test valid saved partitions and round-trip composition. Test high-MOI union deduplication separately. Do not enable unsupported high-MOI NT analyses or preserve malformed internal state solely for a test. |
| I4 — P2 | prepared_index_factory is a current lean lookup closure; native_observed_positions merely casts integers, while actual base conversion occurs in the Rcpp adapter. Legacy calibration/source-cell helpers exist only for test calls or optional worker fallbacks. | Rename the real closure by purpose, remove test-only aliases and reconstruction fallbacks after migrating actual callers, and expose the real one-based R/zero-based native contract. make_fit_view gets exact source/analysis/fit maps, order/duplicate rules, NULL-group semantics and a numerical example. Do not materialize a per-pair/all-guide view bank. |
| E1 — P2 | Drivers take 24/23 formals plus a context containing the same data. Nullable context/index/row/cell fallbacks are unused by actual package/consumer calls. Unread QC thresholds and repeated cache lookups remain. | Apply reviewed E1–E7 interfaces: one ordinary prepared context, small execution options and explicit scientific workers. Require prepared inputs, direct calls, one cache lookup and one common-design/view preparation owner. Preserve permutation response traversal, CRT response prepass/guide traversal, NT group fits, full-B bank scope and cache insertion order. |
| E2 — P2 | get_id_from_idx combines ID lookup, progress/file writes and periodic GC; mixture assignment imports it from association execution. | Plain lookup in callers; separately named reporting and periodic maintenance in execution.R. Preserve exact output bytes, labels, first/every-five reporting, every-200 GC even with progress off, and CRT prepass/guide log append lifecycle. Test observable messages/files and cadence. |
| E3 — P2 | Mixture use_glm=FALSE is uncalled and broken; only TRUE reaches the active EM path. Coefficient-only assignment fits/reduced score mode also have no production caller, only aliases/new tests. Maximum assignment takes an entire S4 object unnecessarily. | Remove unused private modes/flags, pass explicit matrices/metadata to maximum assignment, and name retained fit/score outputs plainly. Retain actual low-count/nonconvergence fallback, warning policy and starting-guess/RNG behavior. Tests do not create a feature obligation for unused modes. |
| D1 — P2 | Old precomputation aliases, construct_data_frame_v2, unused native score/effect wrappers, compute_colwise_max, all check_functions.cpp diagnostics, and production test name setters have no genuine runtime/consumer need. | Delete after per-symbol call readback. Move the used dimname operation into test support; remove the unused one. Do not regenerate frozen inputs. Keep process_calibration_result/process_discovery_result/apply_grouping_to_result and real pipeline helpers as thin, documented adapters. Retire the unused slow CRT sampler from production; preserve its historical oracle/config/archive checks separately from the three active candidate-versus-reference full-bank comparisons. |
| D2 — P2 | Version/history names and the generic shared native header obscure ownership. Only the stored-position loader is shared with count kernels; tail helpers belong to inference. | Produce a retained-symbol name crosswalk, migrate actual callers/bindings atomically, and remove obsolete aliases. Place inference declarations in inference ownership and the sparse loader with count/storage ownership. Public names and genuine component contracts remain unchanged. No blanket rename without call evidence. |
| O1 — P2 | Triplet conversion can reorder caller-owned sparse slots; construct_trans_pairs changes a supplied data.frame to data.table by reference. These are accidental ownership effects, not numerical algorithms. | Proposed bounded ownership corrections, explicitly recorded as observable side-effect fixes: owned sorting vectors/local pair table, identical returned values/classes/order/dimnames. Triplet cost is estimated at roughly 16 transient bytes per stored entry, no dense realization; [the conversion-resource protocol](cleanup/triplet-resource-plan.md) measures total sorting workspace and peaks before acceptance. Verify against independent baseline copies and caller serialization; review the cost rather than add ownership heuristics. |
| C1 — P2 | Most internal functions lack useful contracts; numbering/comments are stale, public plot return prose is sometimes reversed, old class comments/unused locals remain. Even the new cache contract incorrectly labels the stored source-row design “analysis-cell order.” | Review all retained named functions against actual callers. Concise internal source contracts cover shapes/coordinates/order/duplicates/effects/cache/RNG/files and reasons for unusual cases; brief/shared contracts for trivial helpers, useful comments for local workers, no boilerplate for anonymous callbacks. Correct public help text in a narrow separate docs diff with unchanged usages/aliases. AST parity is evidence of unchanged code, not of correct prose. |
| C2 — P3 | Processor validation repeats in three facade paths; metadata-reader connections close only on success; miscellaneous dead locals/comments obscure control flow. | Consolidate the actual shared validation with identical normal/error ordering, guarantee file cleanup, remove proven dead state. Do not use this as blanket permission to redesign all invalid-input diagnostics. |
| B1 — valid-input bug | Unsorted accepted coordinate tables can yield malformed CSC pointers; Cell Ranger modality remapping assumes contiguous feature blocks. Sorting/coalescing choices affect order and duplicate semantics. | Adopt the bounded ordering/modality correction under the ledger below. First characterize grouped/unsorted/duplicate/explicit-zero and interleaved-modality cases. A correction must define stable ordering and duplicate policy explicitly; never silently substitute sparseMatrix coalescing. |
| B2 — valid-input bug | Native cis construction exits after the first in-range block although custom response positions need not be sorted; e.g. 90,1000,110 around 100 omits 110. | Adopt the bounded correction below: remove the early exit and preserve input row order; add the unsorted custom-table regression and safe coordinate arithmetic/preconditions. Pair families/results can change for affected inputs. |
| B3 — valid-input bug | NT-only intercept designs currently fail; drivers deliberately reproduce dimension dropping. Default multi-column fixtures do not justify this as a good design. | Adopt the bounded matrix-dimension correction before or atomically with CL3, with an independent numerical test. Do not implement deliberate vector/error compatibility machinery. See the adopted/deferred ledger below. |
| B4 — valid-input plot bugs | Nondefault QC plot options are ignored; guide bins can take log2 of a negative value; significance-color boundary, Bonferroni guide/target annotation and facet QQ selection have confirmed or source-inferred inconsistencies (L13–L20). | Keep confirmation separate from inference. Scope fixes individually with reproductions and intended visual semantics; do not treat them as documentation-only or regenerate scientific outputs to make plots pass. |
| B5 — invalid-input/output policy | Wrong mixture upper-range variable, validator error precedence, broader missing alpha/QC checks, dispersion catch-all masking nonnumerical failures. These differ from the S1/S2 memory-safety prerequisites. | Explicit dispositions below. Preserve numerical fallback arithmetic; characterize failure classes before narrowing catch policy. Additional public rejection/error changes are not a mechanical rename. |
| B6 — storage/feature asymmetry | Vector-ID in-memory aggregation densifies counts; backed import does not perform the same aggregation. | A sparse replacement may preserve existing in-memory mathematical results, but backed aggregation would change supported behavior. Characterize the actual aggregation/duplicate/output contract and resource profile first; no new backed feature under a cleanup label. |

## 4. Target design and concrete change paths

Use the existing package/modules. No new classes, backend registry, plugin API,
generic test engine, serialized schema or extension hooks are needed.

- The S4 facade validates supported requests and owns commits. Its transient
  context owns canonical maps, requested rows, one common fit view/design and
  the existing compact cache. The saved source design remains a source-row
  matrix; cache validity compares it conservatively. Workers receive aligned
  ordinary inputs and never rebuild an index from raw assignment arguments.
- Sparse/backed operations belong to storage. Owned count kernels consume
  prepared coordinates and own native base conversion; only ondisc adapters
  localize/copy for its current API. Import-format parsers establish their
  representation assumptions explicitly.
- Fitting, score preparation and native inference retain current mathematics.
  Remove unused modes instead of providing speculative fit registries. Keep
  the active mixture fallback and staged-tail behavior as defined algorithms.
- Execution retains two plans and four actual fitting cases. The concrete
  signatures, caller/deletion map and lifecycle risks are in E1–E11. Reporting
  and GC are explicit worker effects, independent of ID lookup.
- Plain result functions retain guide integration, failed-QC rows and complete
  family correction. Real pipeline wrappers stay at that boundary; unused
  forwarding names disappear. Nuclear artifacts remain raw.

Necessary retained complexity must have a named current reason: 56-slot saved
objects; low-MOI NT-local saved memberships; distinct permutation/CRT reuse;
per-guide NT fits; ondisc private interfaces; stage-specific draw segments and
fallbacks; active mixture fallback; actual pipeline part/result contracts.
“The old code did it” and “a test calls it” are not sufficient reasons.

The concrete symbol migrations and retained consumers are enumerated in
[the private-symbol crosswalk](cleanup/symbol-dispositions.md); no aliases are
added to satisfy tests alone. The reviewed target retains the useful lean
treatment-lookup closure as `make_treatment_index_lookup`. Core audit C05’s
earlier factory-removal proposal is superseded. Execution E7’s deliberate
one-column drop/error preservation is also superseded: B3 keeps design matrices
as matrices before or atomically with context migration.

The concrete native ownership target is one sparse-count implementation unit
for the stored-position loader, pair summaries and calibration group summaries.
The loader is private to that unit; it does not belong in an all-purpose
statistical header. Move the group-count kernel out of calibration sampling.
Its Rcpp boundaries take owned coordinate vectors, translate one-based R
positions once, and leave caller vectors unchanged. Delete the superseded
pairwise/shared-count implementations after actual callers migrate. Inference
owns its cross-file tail declarations; tail-local fitting/check helpers stay
private there. Rcpp bindings are regenerated, never hand-maintained.

Native declarations need the same substantive contracts as R functions:
array dimensions/order, coordinate bases, draw-bank bounds/lifetime,
mutation/ownership, fallback semantics and returned records. Document these at
the owning declaration or component, with proportionate local explanation;
do not duplicate ceremonial parameter lists throughout the implementation.

## 5. Behavior scope requiring an explicit disposition

S1/S2/S3 protect native memory/allocation contracts while leaving ordinary valid
draw streams unchanged. O1 removes accidental mutation with the same returned
mathematical results and a measured sparse-copy cost. These are named bounded
corrections, not new statistical functionality.

The senior reviewer recommends the bounded corrections below as engineering
repairs within the existing requested behavior. They select no new statistical
method. The original normal-case fixtures remain unchanged; new regressions
cover the repaired domain independently. This is an explicit exception ledger,
not a blanket claim of identical behavior on every previously accepted input.

| Domain | Disposition and boundary |
| --- | --- |
| B1 imports | Adopt stable cell/feature ordering and actual modality-ID-to-local-row mapping. Preserve stored entries, explicit zeros and the characterized duplicate policy. No coalescing/summation policy is selected. Stop if representation validity cannot be repaired while preserving that policy. |
| B2 cis pairs | Adopt the stated distance predicate for unsorted supplied responses, with safe arithmetic and supplied row order. Affected pair families can change; no new distance method. |
| B3 one-column designs | Adopt `drop=FALSE` for already accepted intercept-only formulas. Test NT calibration/discovery and both traversal modes with independent aligned calculations. Keep saved full-source design/cache authority. |
| B4 plotting | Adopt documented option forwarding. Adopt bin-boundary, existing-significance-flag coloring and one-band-per-panel repairs only after reproductions establish the existing stated contract. Defer a new Bonferroni guide/target interpretation and new global downsampling cap. Retain/document current per-bin sampling, RNG and grouped return class. |
| B5 diagnostics/fallbacks | Adopt wrong mixture-range-variable correction and guaranteed connection cleanup. Defer broad public validation/error precedence changes and changed numerical fallback policy. Characterize actual caught numerical exception types and returned state; preserve their suppression/arithmetic while allowing unrelated allocation/system failures to propagate as a separately reviewed bounded exception-policy correction. Stop if that taxonomy requires a new numerical policy; do not silently change MLE/MM semantics. |
| B6 vector aggregation | Retain the current backed capability boundary. A sparse replacement of existing in-memory aggregation is allowed only with exact order/duplicate/rounding/output and resource evidence. New backed aggregation is deferred new functionality. |
| O1 ownership | Adopt caller nonmutation with identical returned mathematics/schema/order and independent original copies. Measure triplet sorting workspace/peak RSS/timing; review an explicit bounded correctness-cost exception if needed. The 16-bytes/entry estimate is not a measurement or a passing gate. |

The complete per-finding mapping, including L01–L20 and reproduction prerequisites,
is in [behavior-dispositions.md](cleanup/behavior-dispositions.md). Deferred
items are known limitations with named current reasons, not erased findings.
They do not block the agreed internal refactor; a new requirement to change
scientific fallback, guide/target interpretation, duplicate policy or backend
capability would require Eugene's explicit choice. There is no such new method
or feature in this plan.

### 5.1 Exact sampler/calibration safety contracts (CL1)

These checks precede unsafe allocation/indexing. R-facing scalar counts must be
length-one, finite, integral and representable before integer conversion; C++
integer entrypoints reject negative/NA sentinels. Use checked wide sums/products
before narrowing to existing int-indexed storage or R matrix dimensions. Bank
construction holds its heap object with RAII until successful external-pointer
registration, so a later allocation/error cannot leak the partially built bank.
No new pointer class or serialized format is needed. Check
actual container/R allocation capacity, not an arbitrary new cardinality cap.
Retain supported safe zeros. Ordinary valid requests keep operation/draw order,
Boost versus std generator families, fixed seeds and all numerical arithmetic.

| Boundary | Required predicate, placement and disposition |
| --- | --- |
| Calibration request | After assignment/QC and before combination construction: use the actual current nonempty `grna_assignments$indiv_nt_grna_idxs` groups (not guide metadata). Require `1 <= calibration_group_size <= active_nt_count`. Do not clamp the requested/default group size. Require positive integral `n_calibration_pairs` representable by its native consumer before derived B3 conversion. Keep the existing zero-pair diagnostic. Test a whole NT group lost during QC. |
| Derived budgets | Validate finite representable `ceiling(mult_fact*n_pairs/alpha)` before narrowing, and nonnegative B1/B2/B3 with checked sum. These checks protect the existing allocation contract, not a general redesign of all parameter validation. Stage 1 can legitimately use B2/B3 zero. |
| Complement calibration M | Compute the old largest-group sum without integer overflow; use `min(old_upper_bound,n_cells)` only for complement calibration. High-MOI union size cannot exceed this bound. NT-control low-MOI calibration retains its existing M and validated supported partition. |
| Fisher–Yates producer | `n_tot >= 0`, `0 <= M <= n_tot`, `B >= 0`. M=0, B=0 and n_tot=M=0 are safe and remain supported. Validate before vector creation; all formerly safe bank contents are identical. |
| Hybrid producer | Actual algorithm requires positive control and minimum treatment sizes: `N >= 1`, `1 <= m <= M`, `B >= 0`; checked `N+M` and `M+1` must fit the existing native indexing/allocation model. Public NT discovery already requires nonempty controls/groups. Do not invent m=0/N=0 behavior; B=0 remains supported for otherwise valid dimensions. |
| Active CRT producer | `B >= 0`; every fitted probability finite and in [0,1], checked once per constructed bank. Probabilities 0/1, zero-length probability vectors and B=0 are safe. Empty individual CRT draws remain supported. Never validate by redrawing or rescan all bank entries per tested pair. |
| Combination enumeration | Actual current NT count and `1 <= group_size <= n_nt_grnas`; nonnegative integral enumerated count equal to the represented combination count and allocatable matrix shape. The caller enumerates only the existing <=100 branch. Preserve next-permutation ordering. |
| Combination sampling | Same active NT/group bounds; positive finite pair/gene counts and `0 < p_hat <= 1`. Informative error for zero/nonfinite rate: no probability floor/new sampler. Compute the existing ceil estimate in double, clamp against its existing threshold/possible-group bound before checked narrowing; reject an unrepresentable/nonfinite selected count. A +Inf combinatorial upper bound may remain a saturated bound, never an allocation size/cast. Preserve std RNG/unordered-set output order. |
| Candidate-pair sampling | Matching group-count/total gene dimensions; nonempty group definitions with valid zero-based group rows; positive integral gene/request counts and nonnegative count prerequisites. Check `n_groups*n_genes` in a wide type before int allocation/indexing and check count sum/subtraction representability before arithmetic. Retain the full existing candidate permutation, seed, QC predicate, early stopping and one-based returned IDs. No new streaming sampler. |
| Native test adapter/segments | Aligned y/mu/a/w and D columns, observed bounds and `n_trt` agreement; nonnegative budgets/checked segment end against bank length. Bank producers own coordinate validity once; each selected draw checks only its required prefix length (or uses its actual CRT length) in the existing draw loop. No full bank-content scan per pair. Trusted internal components document the producer precondition. |
| Skew-normal fitter | Require at least 100 statistics when the skew-normal fitter is actually reached (largest checked quantile .99); fail before fitting/indexing on a short reached segment. Do not reject unused short B2 metadata after a stage-1 result. Preserve normal supported tail math, sort/reverse behavior and B3=0 fallback distribution. |

Tests cover safe zero banks/samples and p=0/1, oversized/negative/NA requests,
post-QC NT group loss, overlap-heavy complement calibration, zero empirical
candidate rate, representability failures without huge allocations, short
reached skew-normal versus short unused metadata, all frozen full-bank streams,
and all 17 native inference cases. No unsafe baseline case is executed.

## 6. Migration and bounded review units

| Unit | Scope/dependencies | Acceptance before next dependent unit |
| --- | --- | --- |
| CL0 | This inventory/register, source pins, safe WIP checkpoints, actual-consumer list, independent senior review of the plan | Every material item classified; P0/P1 plan gaps resolved; important behavior choices explicit; no further broad edits before review |
| CL1 | S1/S2/S3 sampler/calibration safety and specific native preconditions, before broad restructuring | Safe overlap/invalid-boundary regressions; all ordinary frozen full-bank RNG/scientific cases unchanged; no execution of old unsafe examples |
| CL2 | I1/I2 canonical nonmutating count boundaries; O1 ownership fixes as adopted below | Apply/rework checkpoint atomically with bindings/actual callers/tests; sparse/backed counts, NT sums/order/duplicates, caller nonmutation and bounded memory verified |
| CL3 | I3/I4 + E1/E3 and B3 context/view ownership, required matrix design/driver/worker inputs, real mode pruning, private naming/deletions | Coherent context+worker migration, no temporary alias/fallback layer; API/saved/workflow/cache/native fixtures and fit/draw call scopes pass |
| CL4 | E2 execution effects, remaining dead-code/header ownership, shared validator/file cleanup | Exact messages/file lifecycle/GC cadence; real mixture success/fallback and serial/fork behavior; no production test residue or unowned declarations |
| CL5 | C1/C2 source contracts/comments/public help, reviewed checkpoint adoption; all retained modules | Per-function checklist completed and text checked against callers; source AST/public API checks appropriate to docs-only units; generated help whitelist |
| CL6 | Only explicitly adopted B corrections, each independently reviewable | Reproduction, intended result/error/ownership contract and separate new regression; unchanged normal-case baseline; no broad snapshot replacement |
| CL7 | Integrated source review and qualification after final code | All findings closed or justified with actual boundary/consumer reason; final-source evidence and honest remaining limitations |

Units can split by coherent ownership and run in isolated worktrees. Context
and worker changes may need one atomic commit; do not create fake compatibility
just to split them. Named correctness changes stay distinguishable from
structural/docs changes. Preserve recoverable checkpoints; main is not merged.

## 7. Verification and stopping rules

For each change, test its real contract, not just an implementation mirror.
Keep all frozen source/input/expected artifacts unchanged. Numeric tolerance
stays 1e-8; IDs, order, factors, budgets, QC/significance, stages and supported
errors/RNG are exact. Added bug regressions are separate from frozen results.

- Index/count tests: noncontiguous source cells, ordered NT and group views,
  repeated/unsorted requests, overlapping high-MOI membership, zero/fractional
  stored entries, empty groups, independent reopened backing files, caller
  serialization and the explicit make_fit_view example.
- Workflow/component tests: public/formal/method/slot inventory, saved objects,
  warm/cold cache mutation, full native cases/samplers, real mixture fallbacks,
  serial/two-worker paths and actual fit/bank generation counts.
- Qualify the final immutable install with the verified original in the same
  release toolchain/dependencies; full installed suite, native components and
  paired package checks, plus Linux/macOS/Windows/required-backed/no-ondisc CI.
- Repeat affected actual consumers after code changes: 30 stage scripts and
  full standard/nuclear Nextflow under a shared arrival/work plan. Earlier
  28df390 runs found baseline factor/tie ordering dependent on gather arrival;
  fair=true plus serialized local execution establishes the deterministic
  control. Preserve failed controls and all raw-pod evidence; never sort to hide
  a mismatch. Book coverage remains 14 completed cells plus the same empty-
  barcode prerequisite failure, not a successful whole-book claim.
- Requalify final resource behavior using matched fresh -O2 builds, fixed
  threads/data/budgets/workers, cold/warm repetitions and process-tree RSS.
  Preserve the earlier isolated aggregate outlier/nonreproduction, rather than
  erasing it. Watch owned sparse-copy costs and avoid all-pair view banks/dense
  assays. Static metrics are refreshed after stable cleanup, with costs shown.
- Senior final source review checks actual caller migration, deletion,
  responsibilities, contracts and retained-complexity reasons. Metrics, a guide
  and passing fixtures alone do not establish completion or bug freedom.

Resolve routine failures autonomously. Stop for an unresolved important scope/
scientific behavior decision or unexpected blocker, with the exact conflicting
contracts and a concrete choice. Do not declare completion while S1 or a
material review gap remains.

## 8. Review record

Independent senior review of the first draft: **changes requested**, no P0;
P1 post-QC NT cardinality plus precise guard/crosswalk/scope/resource amendments.
The original signed/hash-pinned review is retained in
[senior-review-initial.md](cleanup/senior-review-initial.md). Sections 4–5 and
the companion crosswalk/disposition/resource plans address those requests.
**Amendment review pending; cleanup implementation remains paused.**
