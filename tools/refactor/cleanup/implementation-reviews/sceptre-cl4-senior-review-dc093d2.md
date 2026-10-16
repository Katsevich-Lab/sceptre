# CL4/C1 independent source review

Exact candidate `dc093d27acefd2b3897eaf8b5a0bb4e9cbe7282f`, base `d3746f6f9aecd0dfb6bf954ce982dd81659b0026` in `refactor-rf2-storage/sceptre`. Read complete retained implementation/header contracts, changed R callers, generated bindings, tests/oracle/ownership probes and source-static tooling. Independently ran the lightweight Python/lizard lexical check only: 49 retained bodies identical after the approved crosswalk, ten exact retirements, signatures identical except the private read-only outlier const qualifier, three R callers changed only by renames/two unread-local removals, protected sources unchanged. No R process, compilation, installed test or benchmark.

**Outcome: no P0/P1 or new material P2 introduced by this rename/retirement/comment unit. Two pre-existing material P2 native lifetime findings remain open below. This source review is not global ownership, ABI/link or runtime qualification.**

## Coherent cleanup

- Actual association/calibration/CRT names replace version suffixes atomically in definitions, the four R workers, calibration constructor/facade, both generated bindings and executable test/component/probe lookups. Retired slow CRT and its generated registration disappear, with no aliases or fake fast replacement. The ten removed bodies match the signed disposition; other real representation primitives remain.
- Tail public-to-module declarations now live in inference.h; local fit/outlier/goodness helpers use an anonymous namespace. Namespace resolution between score/stage/tail/entry sources is coherent on inspection. Numerical expressions, selected-stage offsets, <=.02 transition, B3=0 reuse and returned record schema are unchanged. Compilation still needs to verify linkage.
- Generated R wrappers and .Call registrations agree in name/arity, preserve RNGScope and existing invisibility, and remove old registrations. NAMESPACE/DESCRIPTION/public formals and frozen sources/artifacts are unchanged by this unit. Static hashes support this claim; no regeneration was performed by this reviewer.
- Three active samplers use a literal candidate-name map and installed namespace lookup; all requested rows are materialized, including 5498 fast-CRT rows. Expected values come from the unchanged original-source bridge. Exact four historical configs/source/archive hashes are asserted on every platform. The retired slow sampler has a separate recorded-platform oracle/archive comparison; no symbol-existence filtering weakens active coverage.
- The historical generator now calls verified original names locally, rather than accidentally using candidate helpers. It retains strict original-source verification/clean compile/reverification. No fixture was regenerated. Ownership probe migration correctly changes four old producers to the three retained producers and nine controlled pointer-boundary cases. These cases do not test all native return allocation.

## Maintainer contract review

The comments explain the actual algorithm instead of only naming wrappers. Checked fixed-size extraction order and hybrid m<=k<=M prefix universes; CRT column-wise binomial/subset stream and variable cardinalities/full-bank sharing; R observed versus native synthetic bases; score/effect dimensions/formulas; inclusive/add-one sided tails, unclamped two-sided values, asymmetric skew cap, goodness/outlier cutoffs, p floor and observable sorting/sign reversal; calibration std-versus-Boost engines, combination iteration/unordered-set order and sum/union responsibilities; stored-position counts versus loaded-value >.5 summaries; first-max/zero-total behavior; owned-versus-borrowed columns/materialization/allocation boundaries in the already corrected count scopes. They agree with the actual expressions and known callers.

The scope is 49 handwritten implementations in 13 files, excluding separately reviewed dispersion/EM and generated code. Nine interface declarations/types and the cis lambda are explicitly separate. This is a useful navigation inventory, not a count proving source quality. Existing source-order dependencies and actual per-guide/per-response sharing remain visible.

## Open material findings, predating CL4

1. **P2: calibration result R allocations can bypass native destructors.** `src/negative_control_functions.cpp:74` owns the unordered set and sample vector when `IntegerMatrix out` allocates at98. `sample_calibration_pairs` owns five native vectors at194–195 when `List::create` at240–243 allocates/wraps output vectors. Raw R allocation longjmp can skip their destructors. Renames did not create or resolve this. Use the narrow owned-result lifetime repair now assigned separately, preserving candidate streams/order/schema; controlled raw-R-error and C++-failure proof must include final generated-wrapper boundaries.
2. **P2: inference result wrapping occurs while TailResult owns a native segment.** `src/low_level_full_test.cpp:53–57` retains tail.null_statistics across `association_result_record`; `src/inference_stages.cpp:68–84` constructs R numeric/named-list outputs and optionally wraps the segment. Raw R failure can strand that storage for output1/2 as well as3. Preserve stage/math/sorted segment while repairing the allocation boundary; do not call NULL-XPtr-producer tests proof of this path.

These findings were independently cross-checked with the fresh whole-package reviewer. The already controlled CL2 count scopes and NULL-XPtr-before-bank producers remain bounded valid reviews. Neither implies that all native exceptional lifetimes are safe. In particular a blanket outer Rcpp::unwindProtect callback is not automatically sufficient: callback-local C++ owners can be jumped over, and this installed API allocates its continuation token before setjmp.

## Required integration readback

Root's active `tools/refactor/resource-allocations.R:41` still names the original four samplers at this branch pin and must receive its explicit original/candidate map before final diagnostics; this was predeclared root-owned work, not silently historical evidence. Final ARCHITECTURE/metrics/CI claims must use final names; historical plans/configs/results keep old identities. Full-bank comparisons intentionally require an installed package, which contributor instructions must make clear.

Run the fresh immutable installed native/components/full-bank/index/calibration/context/API/saved-state gates after coherent integration and the separately assigned lifetime repair. Keep B6 precision/resources and the accepted numeric triplet cost disposition distinct. No compiled/runtime gate is inferred from this source-static pass.

## Evidence binding

- `R/RcppExports.R` SHA256 `4c86764ee0dea6e769474ce3a84d18c055ce11fb81980c535b451121091082a7`
- `src/RcppExports.cpp` SHA256 `e03567cbecc5799589ae43000663b287dc37ebcfe0423a8258fd7a593dcfbff3`
- `tools/refactor/cleanup/native-closure-coverage.json` SHA256 `f3b11ae0b2b6134bf696ee18baa151dfd20d041e6e0ddf13f2f76bb6626d8687`
- `tools/refactor/cleanup/native-closure-review.md` SHA256 `a3e66a0c616679801404cc2d2ad53a1fa2abaff229fc19af769bb336bafc3b3e`
- `tools/refactor/native-closure-static.py` SHA256 `7bb742259102e212e4d6816428d6d81d410c857346bdda8b2f4a4aa66a660e3a`
