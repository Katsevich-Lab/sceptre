# Refactor status — October 16, 2026

RF0 is the frozen baseline and characterization checkpoint. Implementation is
under way on `codex/refactor`, which is published for hosted CI. The accepted
scope remains identical sceptre API and no new functionality. Nothing has been
merged to main or released.

## Implementation progress

- RF1: workflow prerequisites/transitions have one owner, retaining the saved
  flags and sibling analyses. Deeper 12-cell source/analysis/fit/native boundary
  tests passed against the untouched reference before migration. Index/context
  implementation is in an isolated worker worktree.
- RF2/RF3: the selected route is a private sceptre storage adapter, with no new
  ondisc API or release prerequisite. Storage extraction is in a second worktree.
- RF5: current response-cache invalidation and facade writeback have a single
  module owner. Actual design matrices (including contrasts) and control views
  govern validity. PR 233 remains open with changes requested at its unchanged
  head `3f1ff0650078b382662ee70b72ea5c72fa15b031` (verified October 16); its new
  guide cache is deferred and is not being recreated here.
- RF7: native inference extraction is in a third isolated worktree.
- The full refactor goal is active. Routine implementation/test failures are
  resolved autonomously; important unresolved design choices or unexpected
  blockers are the stopping conditions. Main merge/release remain separate.

### Bounded cache correctness correction

An untouched-baseline experiment confirmed that NT-only calibration followed by
changing to complement controls retains four incompatible response fits and
produces warm/cold result differences. The cache-validity contract in section 9
requires clearing those fits. The refactor now invalidates on control-view or
actual-design changes, without changing the API or adding cache capabilities.
Regression tests check that correction, contrast changes with identical formula
text, and discovery→calibration→discovery warm/cold parity in three workflows.
The frozen supported-workflow, cache, saved-object and API fixtures still pass.

## Source and environment

- sceptre source: `21f9ea098b69c444a884a49bc254de711968e3b5`, 0.99.0,
  freshly fetched `origin/main`; PR 233 excluded.
- Selected ondisc source: 1.3.5, upstream `589803edea2a66ad15e7f07186bee16c7739ab8c`,
  same source tree as local `4d9d75e18f21d11b6da6a54cb03d4424fcf16fe4`.
  Built from a Git archive into a temporary library; no user-library replacement.
- Earlier environment observation: installed ondisc 1.2.0. It produced four
  bounds warnings in the synthetic backed preparation. The 1.3.5 source lane
  removes those warnings. These are separate dependency observations.
- Pipeline and book revisions, actual scripts/slots/containers, observed runtime
  availability and existing limitations are in `consumer-inventory.json`.
- R 4.6.0 on macOS arm64; full dependency/platform/BLAS/RNG details and fixture
  hashes are in `tests/testthat/fixtures/refactor/manifest.json`.

## Evidence established

- Before RF0 test additions, the unchanged package passed 2,578 assertions in
  56 test blocks with no failures, warnings or skips (installed ondisc 1.2.0).
- At RF0, source checks showed no changes to production `R/`, `src/`, `DESCRIPTION` or
  `NAMESPACE`. The S4 facade retains its existing 56 slots.
- Versioned public API, workflow, ordered-index, cache, saved-object, native and
  sparse/backed fixtures are generated from verified baseline source. The
  generators force clean native compilation and recheck source hashes.
- Native characterization includes 17 deliberately reached inference cases and
  four sampler banks, including all 5,498 draws in a production-budget CRT bank.
- The source audit inventories all 11 tracked pipeline R scripts and their
  explicit component dependencies. No claim of complete project knowledge is
  based on this bounded inventory.
- The original source package's real scripts passed 12 fresh-process stage
  smoke tests with ondisc 1.2.0. Two repetitions produced exactly equal
  comparison payloads. This includes standard calibration/power/discovery
  workers and gathers, nuclear QC and trial-mode raw Parquet shards.
- Same-environment original-source `R CMD check --no-manual --no-vignettes`
  completed with two pre-existing vignette warnings because the source build
  deliberately omitted built vignette output. This is not a clean full-vignette
  check; the candidate must use the same flags for comparison.

## Candidate checkpoint verification

- Full suite with required source ondisc 1.3.5: **3,265 assertions in 72 test
  blocks; zero failures, errors, warnings or skips**. This includes the exact
  all-platform sampler oracle, partial-guide QC, real assignment mutation,
  serial/two-worker equivalence and backed reattachment/mismatch checks.
- All 23 recorded fixture SHA256 digests were verified after final generation;
  package tests also verify the saved MD5 inventory.
- Reference and candidate each passed all 12 fresh-process pipeline stages
  with source ondisc 1.3.5. Their complete comparison payloads were **exactly
  equal** (`all.equal(..., tolerance=0)`). Both preparation reports had no
  bounds warnings. Reports/payloads are retained under `results/`.
- Reference and candidate source checks with ondisc 1.3.5 both completed with
  exactly the same two omitted-vignette-output warnings and no errors or notes.
  The warning blocks were compared directly; logs and the comparison are saved
  in `results/package-check-*`.
- The resource harness executed 24 fresh serial R processes: warm-up plus three
  measured repetitions for cold/warm sparse permutation, CRT and backed
  permutation. Median analysis times were 0.083–0.520 seconds and observed
  serial-process peaks were about 267–316 MiB. These small timings establish a
  runnable baseline harness, not acceptance of realistic performance gates.

See `results/tests-summary.json`, `results/benchmark-reference.json` and the
pipeline reports for machine-readable RF0 evidence. Hosted CI run
[38917674952](https://github.com/Katsevich-Lab/sceptre/actions/runs/38917674952)
is running against the published RF0 checkpoint; no hosted success is claimed yet.

The maintainability baseline is frozen in `MAINTAINABILITY.md` and
`results/maintainability-baseline/`: fixed change-tracing exercises and a seven-
operation testability inventory, reviewed ownership sites, plus versioned R/C++
complexity/nesting and candidate-site inventories. The unchanged candidate has
identical detailed metrics and source hashes. Human comprehension timing,
cognitive complexity, token clone coverage and a complete dynamic dependency
graph are explicitly unmeasured; there is no composite or monotonicity target.

## Independent review dispositions

Final independent senior review of the frozen checkpoint found no remaining
P0, P1 or material P2 issue for RF0. The reviewed tree was preserved exactly
while the final unpushed coordinator batch was split into workflow fixtures,
execution evidence, CI and maintainability commits. This review does not mark
the explicitly outstanding hosted/runtime/resource gates complete.

- A generator could otherwise reuse a newer stale native library despite
  matching source files. Both generators now force rebuilding and verify source
  hashes again after compilation.
- The assignment-mutation input now has counts between thresholds 40 and 45,
  changes real memberships and compares the recomputed result with cold work.
- The Bonferroni input now has a target/response with only one passing guide,
  and the test explicitly asserts that mixed outcome.
- The complete fixture manifest must be finalized after native regeneration;
  package tests verify all saved artifact digests.
- Native sampler parity now compares every full bank with a frozen baseline
  source oracle compiled in the same test environment. It runs on all platforms;
  the environment-qualified stored stream is an additional check. A separate
  deliberate draw-order mutation was detected by the mandatory oracle.

## Gates that remain distinct

- Hosted Linux/macOS/Windows CI and the explicit ondisc-free lane are running;
  their final results remain a gate. The backed lane must fail when ondisc is
  missing, not silently skip.
- Full Nextflow/nf-test and the historical 0.10.2 containers have not run. The
  current profiles contain NULL placeholder inputs, local Docker has no running
  daemon and nf-test is absent.
- Nuclear nontrial and zero-passing-pair script branches, upstream parameter/
  assignment/pod-planning scripts and an uncached book execution need additional
  consumer coverage before their structural migration. Current smoke scope is
  explicitly trial-mode and fixed supplied pod maps.
- Stored sparse zeros/fractional counts, reordered caller cell requests,
  discovery→calibration→discovery cache reuse and deeper count/score/native-base
  indexing coverage need characterization before their affected RF1/RF3/RF5 moves.
- Small serial synthetic timings are diagnostic. Realistic difficult-fit data,
  a backed workload larger than an explicit dense-memory budget and fork-worker/
  process-tree peak RSS remain unmeasured acceptance gates.
- API/scientific/default changes, mixture repair, replacement QC, accelerated
  fits, SPA, Wilcoxon, backed CRT and new storage backends remain separate work.
