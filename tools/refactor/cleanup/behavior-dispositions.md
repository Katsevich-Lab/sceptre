# Cleanup behavior dispositions

Status: **planning decisions; implementation and qualification pending**.
These decisions amend the open dispositions in [CLEANUP_PLAN.md](../CLEANUP_PLAN.md).
They do not modify the original [L01–L20 findings](legacy-findings.json), their
source hashes, or the [native findings](native-review.md).

Authority: independent cleanup-plan review of source
`c7d62502757f2a17741e2ef18b9e4d006bf063dc`, October 16, 2026. Its original
artifact is `/private/tmp/sceptre-cleanup-plan-senior-review.md`, SHA256
`6262bc716f01070fdfdaab4db267b5c56aa901e048668b77e8476e3c1b390e5e`.
The unchanged legacy register has SHA256
`c8a64b5d8d58ea0cfd547fac085cd158485ac661fb5f6102c45db22db8595a8e`.
Line references in that register remain tied to its original source pin;
function names below identify the affected current boundaries.

“Adopted” means included in the reviewed cleanup plan, not already fixed.
“Reproduction prerequisite” means the bounded correction may proceed only after
the stated reproducer and intended existing contract are established.
“Deferred” names retained behavior and its reason; it does not mean the finding
is absent. Routine adopted engineering repairs require no additional user
confirmation. The explicit stop conditions below cover genuine unresolved
scientific, duplicate-policy, or backend-capability choices.

## B1–B6 decisions

| ID and exact source finding | Disposition, affected domain and reason | Required evidence and constraints |
| --- | --- | --- |
| **B1 — coordinate imports.** `create_csc_matrix()` assumes grouped triplets; Cell Ranger import subtracts a modality's minimum feature ID, assuming contiguous feature blocks. | **Adopted:** stable cell/feature ordering on owned scratch vectors and actual global-feature-ID to local-row matching. Repairs unsorted coordinates and interleaved modality IDs without selecting a new count definition. | T1. Preserve entry multiplicity, tie order, explicit zeros, values, dimensions and IDs; do not coalesce duplicates or substitute a constructor that silently sums them. Preserve the separately characterized duplicate policy. If representation validity requires a new duplicate/rejection policy, stop for that explicit decision. No full dense count realization. |
| **B2 — cis predicate.** `compute_genes_within_distance()` exits after the first in-range block, omitting later matches when supplied positions are unsorted. | **Adopted:** remove the early exit and implement the existing inclusive distance predicate with safe coordinate arithmetic, preserving supplied response-row order. | T2. The supplied positions `90, 1000, 110`, midpoint `100`, distance `20` must select the first and third rows. Check endpoint/arithmetic boundaries before unsafe narrowing. No sort of caller rows, changed distance definition, or new coordinate-NA policy. Affected pair families can change. |
| **B3 — one-column designs.** NT calibration/discovery subsetting drops a matrix dimension although intercept-only formulas are already accepted. | **Adopted:** preserve dimensions with `drop = FALSE` at the real common/group fit-view boundaries, before or atomically with CL3. | T3. Independent aligned calculations for NT calibration and discovery, both permutation and CRT. Keep the saved full-source design and cache-invalidation authority. This supersedes the old execution E7 proposal to encode/preserve the vector/error behavior; do not build temporary compatibility machinery for that error. No new fitting fallback. |
| **B4 — plotting.** QC options are ignored; histogram expansion, volcano colors, Bonferroni IDs and facet QQ bands have separate inconsistencies. | **Adopted:** forward documented QC plotting options. **Reproduction prerequisite:** histogram expansion only beyond `max_single_bin`; volcano colors follow existing `significant` flags; QQ bands select one group per panel. **Deferred:** new Bonferroni guide/target interpretation and a new global downsampling cap. | T4/T5. Confirm each bounded repair independently; do not introduce a new significance/tie rule. Preserve per-bin sampling, RNG, grouped return classes, defaults and scientific results. Bonferroni pooling/labeling is a distinct semantic decision, not a routine forwarding correction. |
| **B5 — diagnostics and fallback policy.** Wrong mixture-range variable, validator precedence/gaps, success-only connection closure, and catch-all dispersion exceptions differ in scope. | **Adopted:** compare the perturbation upper bound to `g_pert_guess_range`; guarantee owned Matrix Market connection cleanup. **Conditionally adopted in this cleanup:** replace the blanket dispersion catch after actual exception/state characterization, preserving numerical suppression while allowing unrelated allocation/system exceptions to propagate. **Deferred:** broad public validation/error precedence and a new numerical fallback policy. | T6/T7/T12. Enumerate numerical exception types from executed BH/Rcpp paths and preserve each existing returned estimate/method state, MLE/MM/pilot arithmetic, iteration order, clipping and warnings. Do not infer the taxonomy from passing normal fixtures. If numerical and unrelated failures cannot be separated without a policy choice, record the concrete conflict and stop that correction. No invented estimator, candidate-rate floor, fallback probability, or method-label correction. S1/S2 safety checks below are adopted separately. |
| **B6 — vector aggregation.** `collapse_grna_matrix()` materializes the R input densely; the R-backed import path does not provide the same vector aggregation. | **Conditionally adopted in this cleanup:** characterize/profile the existing **in-memory** operation, then replace full-input densification with bounded sparse aggregation after T11 establishes contract/resource parity. **Deferred capability:** backed aggregation. | Preserve the current returned dense matrix class, first-seen vector/cell order, metadata-duplicate multiplicity and numerical behavior. Removing avoidable full-assay densification is engineering work in scope; no present profile or speed/RSS improvement is claimed. A conflicting duplicate/rounding contract requires a concrete decision, not silent normalization. Sparse output/backed capability is not implicitly authorized. |

## B7 — subsequently reproduced dense-input conversion loss

Independent source review of da6a079 reproduced an existing P2 in both the
original 21f9ea0 package and the candidate: a valid same-label symmetric dense
matrix loses its lower triangle when automatic structured triplet storage is
copied into a general sparse result. A tolerance-near-symmetric input also loses
its supplied lower value. The safe 2-by-2 direct/public-import reproductions and
source-pinned alternatives are retained in
[the separate disposition](implementation-reviews/sceptre-symmetric-dense-disposition.md).
This addition does not rewrite the original audit or its signed plan artifact.

**Adopted bounded engineering correction:** force general shape through public
Matrix coercion before sparsification, preserving every supplied dense value
and dimname. No new method, missing-value definition, class registry or backed
capability is introduced. Do not use deprecated direct dgT coercion, private
Matrix API or invalid hand-built slots to evade allocation. Sparse/backed inputs
must never densify. Qualify symmetric/near-symmetric direct and public paths
against independent intended values; ordinary rectangular/triangular/integer
cases retain exact output and caller nonmutation. Measure the supported public
route's additional dense-input-only payload copy on matched immutable builds;
retain all attempts and distinguish an explicit reviewed correctness cost from
a general resource-gate pass. Triplet and integrated resource gates stay separate.

## Complete L01–L20 ledger

The source findings remain as recorded in `legacy-findings.json`. This table
disposes of each item, including findings containing both cleanup and behavior
changes. Test identifiers refer to the independent evidence below.

| Finding | Disposition, domain and reason | Test/constraint |
| --- | --- | --- |
| **L01 — missing/partial internal contracts** | **Adopted:** review the paused 39-function draft, complete contracts for every retained named helper, and improve the existing QQ contracts. Remove obsolete functions rather than documenting them into permanence. | T9. Types, shapes, coordinate views, order/duplicates, returns, effects and genuine compatibility reasons; brief/shared contracts for trivial helpers and local workers, no anonymous-callback boilerplate. AST parity does not establish prose accuracy. |
| **L02 — reversed plot return prose/obsolete names** | **Adopted documentation correction:** source/Rd text must describe actual list/combined/single-panel behavior and current assignment names. | T9. Narrow generated-help whitelist; exported signatures, aliases, usages and actual return behavior unchanged. Include the no-discovery-pair single-cellwise-plot branch. |
| **L03 — dead comments/metadata duplication** | **Adopted:** remove commented-out class-union remnants, stale numbering and redundant declarations where no runtime contract changes. Retain explanations of non-obvious algorithms and compatibility. | T9. Preserve Matrix imports, native registration and the complete 56-slot saved schema; apparently unused saved slots are not ordinary dead locals. |
| **L04 — unused locals/formals/repeated guards/theme** | **Adopted engineering cleanup after caller crosswalk:** remove proven dead locals and redundant work, prune unused private arguments with their actual callers, and share identical validation only when evaluation/error order remains the same. | T9 and relevant current consumer checks. No aliases to satisfy tests; no broad validator precedence change. Public formals stay fixed. A private symbol with a real executable consumer retains a named boundary until coordinated migration. |
| **L05 — avoidable forwarding/repeated preparation/version names** | **Adopted:** inline/remove redundant internal forwarding and repeated directory preparation, and rename historical internal names through an explicit retained/deleted-symbol map. | T9/T7 and backed/pipeline artifacts. Preserve actual pipeline calibration/QC parts, restricted/trial guide selection, pod/file helpers and three real result adapters. Preserve native entrypoints until their complete migration; no speculative compatibility layer. |
| **L06 — caller-owned positive-control table mutation** | **Adopted O1 ownership correction:** make local ownership explicit before data.table mutation; accidental caller class mutation is removed while returned pairs stay the same. | T8. Independent original input copies; serialized caller/class/attributes unchanged; returned class, ordering, IDs and values preserved in each exclusion mode. This is an observable side-effect fix, not AST-only documentation. |
| **L07 — dense/backed vector aggregation asymmetry** | **Conditionally adopted in-memory correction:** use B6/T11 to remove full-input densification after exact characterization and resource proof in this cleanup. **Deferred capability:** new backed support. | T11. Preserve sum/duplicate policy, output class and the backed capability boundary. If characterization exposes an actual policy conflict, stop with that concrete choice; absence of evidence today is not a reason to exclude this engineering task. |
| **L08 — scratch CSC and modality assumptions** | **Adopted bounded import correction:** use B1. Sorting and actual ID mapping are explicit behavior repairs, not disguised comment changes. | T1/T8; stable owned ordering, stored entries/zeros/ties preserved, no coalescing. Stop if duplicate-policy preservation and representation validity conflict. |
| **L09 — metadata parser/connection lifetime** | **Adopted resource-lifetime correction:** always close the helper-owned connection. **Deferred:** redesigning accepted malformed-file spacing/header/error behavior. | T7. Existing plain/gzip results unchanged; inject a controlled parse/read failure to check closure, without expanding the file-format contract. |
| **L10 — calibration/QC/alpha validation gaps** | **Adopted S1/S2 subset:** validate actual post-QC active NT group bounds, positive integral calibration requests and representable derived budgets before unsafe conversion/allocation. **Deferred:** general alpha, mito, additional-cell, scalar/NA and unrelated public validation expansion. | T10. Do not clamp group size, count metadata rows as active groups, reject safe zero-bank domains globally, or add an arbitrary cardinality cap. Unrelated invalid-input error behavior remains unchanged. |
| **L11 — wrong mixture upper-range variable** | **Adopted named invalid-input correction:** compare against the perturbation range, as B5 specifies. The current advertised upper bound is the authority. | T6. Boundary/above-bound values with a valid pi range; ordinary mixture assignment/initialization/RNG/fallback paths unchanged. No broad hyperparameter tightening. |
| **L12 — required-column guard/error precedence** | **Deferred:** reordering NA/column validation changes the public malformed-input error contract and is unnecessary for the adopted safety repairs. | Preserve existing error behavior in this cleanup; a future fix needs missing-column/NA/valid-table regressions and an explicit error-precedence decision. |
| **L13 — ignored QC plot options** | **Adopted documented forwarding:** pass the supplied downsample/size/transparency options to the pairwise plot. | T4. Test nondefaults independently and retain default output; scientific QC thresholds/counts/results are not recomputed. |
| **L14 — volcano significance boundary** | **Reproduction prerequisite:** establish the mismatch, then color by the existing `significant` flags, not a newly selected `<`/`<=` threshold convention. | T5. Include the maximum significant p-value and tied p-values; preserve all scientific flags/p-values. Missing/ambiguous flag provenance is a stop for that correction, not permission to invent significance. |
| **L15 — guide histogram expansion bound** | **Reproduction prerequisite:** show a valid high explicit threshold failing when observed counts exceed 10 but not `max_single_bin`; then gate expansion on the actual single-bin bound. | T5. Default/threshold-boundary cases and a high-threshold regression; no alternative binning scheme, new count clipping, or changed assignment threshold. |
| **L16 — per-bin sampling/grouped return** | **Adopted contract correction; deferred global-cap change:** describe per-bin/grid sampling and grouped result class honestly. A global cap would change sampling/output behavior. | T9. Preserve fixed seed, row/bin order, grouping and output class; do not silently ungroup or reinterpret `downsample_pairs`. |
| **L17 — Bonferroni guide/target plotting** | **Deferred semantic correction:** pooling/labels/annotation for guide-level assignments versus target-level results need a defined interpretation. It is not required to simplify the core. | Retain/document the limitation. If a current deliverable requires changing it, stop with a concrete guide/target example and obtain Eugene's choice before implementation. No new interpretation in this plan. |
| **L18 — getter scratch/native ownership** | **Adopted ownership audit/correction where demonstrated:** establish whether the one-guide scratch aliases caller state before pruning its mutating-native route. Prefer an owned/nonmutating boundary. | T8. Raw assignment serialization and repeated one-guide/multiguide getters, both source/QC views, sparse class/IDs/order. No mutation bug is asserted without that evidence. |
| **L19 — production test utilities/mock-validation typo** | **Adopted test-support cleanup:** relocate the actually used dimname helper and retire the unreferenced helper after the private-symbol crosswalk. **Deferred incidental mock-validation behavior change:** it is unrelated to runtime cleanup. | T9 plus existing mock-data consumers. Do not change distributions/seeds, committed fixture inputs, expectations or generation authority; never add a runtime alias solely for tests. |
| **L20 — QQ band global versus per-panel group** | **Reproduction prerequisite:** establish a missing panel band with distinct panel group IDs, then choose one group within each panel under the existing contract. | T5. Multiple panels plus the current single-panel/equal-group-size behavior; preserve rank/tie/scale/band formulas and RNG. No new confidence-band method. |

## Independent regression and resource evidence

Existing ordinary scientific/API/index/saved-object/cache/native/consumer
fixtures remain immutable and retain their original authorities. An adopted
repaired domain needs an independent intended calculation, not expected output
generated from its new helper. If a historical characterization records the
specific defect being repaired (notably the NT intercept-only error), preserve
the original artifact and explicitly supersede only that named case's test
authority. Keep all unrelated cases/gates. Do not broadly skip or replace an
oracle because the candidate differs.

| Test | Required evidence |
| --- | --- |
| **T1 — import coordinates** | Hand-specified unsorted cell/feature triplets and interleaved modality IDs, including explicit zeros and stable duplicate ties. Compare stored entries, pointer/index shape, values, dimensions/dimnames and input serialization against independent ordered expectations. Sorted ordinary files remain identical. Characterize existing safe duplicate cases without using a coalescing constructor as the oracle. |
| **T2 — cis distance** | Compare the unsorted example and inclusive endpoints against `which(abs(as.double(position) - as.double(midpoint)) <= distance)` on supplied rows; include safe wide-coordinate arithmetic. Test public cis construction as well as its native predicate, with supplied row order intact. No unsafe overflow example is run through the old kernel. |
| **T3 — NT intercept design** | Valid, nondegenerate intercept-only data for calibration/discovery and both permutation/CRT traversals. Independently build the source/analysis/fit row mapping and a one-column matrix, then compare fitted quantities/statistics/results. Do not obtain expected rows/design from the new context helper. Existing multicolumn cases and cache scope remain unchanged. |
| **T4 — QC forwarding** | Probe the actual forwarded nondefault options and resulting pairwise layer data/parameters; compare the default branch too. Keep public return shape and existing no-discovery behavior. |
| **T5 — reproduction-gated plots** | Preserve a failing minimal case before each fix: high-threshold histogram; existing-significant flag/threshold ties; two facet panels with distinct group IDs. Assert the intended current contract independently, including single-panel/default behavior, rather than freezing a new whole-plot snapshot. |
| **T6 — mixture bound** | Valid pi range with perturbation upper values at/below/above the documented bound; demonstrate the old wrong-variable path. Existing genuine EM success, low-count threshold fallback, nonconvergence, warnings and RNG cases retain their behavior. |
| **T7 — owned connection/files** | Successful plain/gzip metadata reads plus a controlled failure after opening. Confirm the helper-owned connection closes and returned values/success-path files are unchanged. Do not close caller-owned or unrelated connections. |
| **T8 — ownership** | Run original and candidate on separate independently copied inputs; compare returned mathematics/schema/order/dimnames and caller serialization/class. Include triplet slots, stored zeros/fractions, sparse/backed group counts and repeated assignment getters. Measure the triplet copy **and sorting workspace**: roughly 16 bytes/nnz is an estimate, not measured RSS. |
| **T9 — retained source/contracts** | Per-retained-function caller/contract checklist, complete symbol crosswalk, AST identity for documentation-only changes, narrow Rd text whitelist, exact public API/56-slot schema, and real standard/nuclear pipeline/saved-object/cache behavior. Source text must explain necessary algorithms/compatibility, not preserve dead code because an old test mentions it. |
| **T10 — real safety domains** | Public post-QC loss of a whole NT group, overlap-heavy complement calibration, checked budget/cardinality limits, and native entrypoint-specific boundary cases. Active groups are the current nonempty assignment list. Preserve safe Fisher–Yates zero samples/banks, CRT empty draws and p=0/1, and stage-1 results with unused short B2. Check the skew-normal minimum only when fitting is reached. Do not execute unsafe old combinations/indexing. |
| **T11 — vector aggregation gate** | The bounded evidence specified below, before any sparse replacement. Profiling alone authorizes no output/backend change. |
| **T12 — dispersion exception policy** | Execute the actual BH/Rcpp numerical failure paths; enumerate concrete thrown types and the exact estimate/method state returned at each throw point, including partially completed MLE/MM attempts. Then implement/test a bounded catch policy that preserves suppression and arithmetic for those numerical failures while unrelated allocation/system failures propagate. Controlled test-only fault injection may verify propagation without causing real OOM, but cannot substitute for the executed numerical taxonomy. If the taxonomy cannot be established without a new failure policy, record the concrete unresolved choice. No new numerical fallback or relabeling of historical method results. |

### B6 profiling decision

The in-memory sparse replacement is **conditionally adopted in this cleanup**.
Characterization and profiling are its implementation gates, not a rationale
for leaving full-input densification outside the engineering scope. No profile
or successful replacement is claimed yet. Establish all of the following:

1. Record the actual in-memory `vector_id` contract using independent original
   inputs: dense and supported sparse representations; first-seen vectors;
   original cell/feature IDs/order; single and multiple guides per vector;
   repeated metadata/guide selections; integer, fractional and explicit-zero
   stored values. Preserve the current dense matrix output class and dimnames.
2. Match the original per-vector summation/rounding behavior and duplicate
   multiplicity. Use the existing declared numerical tolerance where applicable
   and exact structure/order/decisions; do not assume reordered sparse algebra
   is equivalent solely because the symbolic sum is the same. A conflicting
   duplicate/sum policy is a stop for an explicit decision.
3. Profile representative, bounded guide-by-cell/vector multiplicities, with
   controlled dimensions/nnz/density and a real input when available. Measure
   allocations, conversion/sorting workspace, peak worker RSS and elapsed time,
   including the unchanged dense output footprint. Avoid deliberate OOM and
   do not present extrapolation as a measurement.
4. Use pinned release builds/dependencies, fixed threads/workers/data, a warm-up
   and at least three measured repetitions. Preserve outliers. Compare median
   elapsed and peak RSS to the original; the existing general gates remain
   median elapsed increase <=10% and RSS increase <=max(10%, 64 MiB). A bounded
   correctness-copy exception must be recorded/reviewed explicitly, never
   relabeled as a pass from an older report.
5. Demonstrate removal of full-input densification with reviewable bounded
   workspace and unchanged outputs, and meet the measured resource gates before
   qualifying the replacement. Investigate a resource/parity failure; do not
   declare the task complete from a profile alone. If exact characterization
   instead exposes a new duplicate, rounding or scientific-policy choice, stop
   with the concrete conflict for review. No backed import/aggregation,
   sparse-output API, streaming backend or new count semantics follows from
   this in-memory correction.

### Safety, retained algorithms and completion limits

CL1's entrypoint-specific table in `CLEANUP_PLAN.md` governs S1/S2. Bounds are
checked before allocation/narrowing against actual active groups and container
capacity, with wide sums/products and no arbitrary cap. A zero/nonfinite
candidate rate gets an informative failure, not a new probability floor or
sampler. Only complement calibration uses
`min(old_largest_group_upper_bound, n_cells)`; formerly safe streams stay
identical. Group size is never silently clamped.

Retain the active Poisson/dispersion arithmetic, mixture initialization and
threshold fallbacks, distinct permutation/CRT sharing scopes, guide-specific
NT fits, compact caches, stored-position versus positive-value count meanings,
backed integer coercion/file identities, and real pipeline/saved-object
adapters. These have current behavioral/consumer reasons; document those
reasons rather than describing every retained algorithm as historical debt.

The adopted corrections intentionally change the named affected domains. This
plan cannot claim universal old behavior, absence of all bugs, or newly
qualified source. Completion requires implementation review and the existing
scientific, exact native-bank, saved/cache, serial/fork, actual-consumer and
affected resource gates on the final immutable build. Earlier 3,522/112/597
results and resource reports remain tied to their recorded builds.

Completion reporting must name the deferred malformed-input precedence/broad
validation, new numerical-fallback policy, Bonferroni plotting interpretation,
global-cap change and backed aggregation capability. It must also report the
actual T12 exception-taxonomy/catch outcome and T11 sparse-replacement evidence,
or the exact unresolved choice preventing either conditional task. An unmet
condition is not a successful implementation. The explicit feature/policy
deferrals do not block the agreed internal cleanup. If
finishing a required deliverable would instead demand a new duplicate policy,
scientific fallback, guide/target interpretation or backend capability, stop
that dependent work with the exact conflicting contracts and a concrete choice
for Eugene; continue independent adopted work. No such new method or feature
is selected here.

## B8–B9 — fresh whole-package review findings

Independent full-source review of cc016835fddaa5dec9a87c9bead61914153433b0 is
preserved in [the report](implementation-reviews/sceptre-whole-package-review-cc01683.md)
and its per-file/flow coverage manifest. These additions do not rewrite the
original signed cleanup plan or claim that its earlier limited reviews covered
these findings.

**B8 — native return lifetime (material P2, inherited).** Calibration group/pair
containers and inference's selected null vector remain live across raw R return
allocation/preservation. The installed Rcpp wrapping path does not automatically
protect C++ destructors from R longjmp. Adopt a narrow lifetime repair with
independent design/source review and safe controlled candidate-source raw-R-
error/native-allocation/success tests plus an intentionally unsafe negative
control. Preserve exact loop/RNG/draw order, return schema, selected/reversed
null segment, reached-stage budgets and exception propagation. Do not introduce
large eager pair/tail capacity merely to move allocation earlier. No real OOM,
catch-all numerical fallback or general ownership framework is authorized here.
Existing CL2/native-bank closure remains scoped to its specifically tested
boundaries; it does not imply this finding already resolved.

**B9 — discovery plotting optional missing diagnostics (material P2, inherited).**
The actual original and cc01683 plotting source both show two valid result rows
and one significant call at output amount 1, but zero rows/calls at amounts 2/3,
because whole-table `na.omit()` requires optional tail/padding fields to exist.
Adopt the bounded consumer correction: determine row completeness from the
actual required plotting fields (`p_value`, `log_2_fold_change`, `significant`),
then retain the same valid scientific family despite missing optional fields.
Missing primary plotting values remain excluded. Stored results, significance,
family correction, sampling/RNG and public formals/defaults remain unchanged.
Independently compare QQ/volcano/text-family data across amounts 1/2/3, including
stage-1/empirical-tail missingness and null padding; verify caller nonmutation.
Do not regenerate scientific fixtures or add a new significance definition.

The same review's P3 changes are ordinary bounded cleanup: remove unconsumed
transient `control_fit_pos`, describe D rows as eigencoordinates rather than
original coefficients, and qualify actual backed nuclear/power guide behavior.
No saved-schema field or future feature follows from those corrections.
