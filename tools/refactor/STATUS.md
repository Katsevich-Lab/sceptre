# Refactor status — October 16, 2026

RF0 is the frozen baseline and characterization checkpoint. Implementation is
under way on `codex/refactor`, which is published for hosted CI. The accepted
scope remains identical sceptre API and no new functionality. Nothing has been
merged to main or released.

## Implementation progress

- RF1: canonical source/analysis/fit maps and request-row lookup are integrated
  into assignment/QC, calibration and both execution plans. The 12-cell numerical
  fixture passed against the untouched reference before migration. Workflow
  prerequisites have one owner, retaining saved flags and sibling analyses.
- RF2/RF3: the private storage adapter is integrated; all six private ondisc
  calls are in `R/storage.R`. No companion API/release is required. Tests cover
  repeats/order, caller nonmutation, stored zeros, fractional coercion, shared
  import identities and reopened files.
- RF4: measurements and policy are separate. Vector eligibility and calibration
  acceptance use the compiled legacy policy; its missing-value compatibility is
  also owned by that header. No QC method or threshold has changed.
- RF5: current response-cache invalidation and facade writeback have a single
  module owner. Actual design matrices (including contrasts) and control views
  govern validity. PR 233 remains open with changes requested at its unchanged
  head `3f1ff0650078b382662ee70b72ea5c72fa15b031` (verified October 16); its new
  guide cache is deferred and is not being recreated here.
- RF5 also separates current numerical fitting and score preparation into
  ordinary modules, with compatibility delegates for the old private names.
- RF6: shared dispatch/logging/binding mechanics are integrated. Permutations
  retain response traversal; CRT retains response precomputation and guide
  traversal, assignment fit and full-bank reuse.
- RF7: native effect/score, staged-tail and result-record components are
  integrated. R observed positions convert once at compatibility boundaries;
  new native components and synthetic banks use zero-based positions.
- RF8: plain table conversion, QC-failure assembly, guide integration and family
  correction have one owner. The facade prepares context once and commits
  returned caches/results. Pinned component names and artifact contracts remain.
- RF9: the superseded index factory was deleted, association modules now have
  responsibility-based filenames, and `ARCHITECTURE.md` provides concrete change
  paths. Final evidence/maintainability reporting and qualification are in progress.
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

## Integrated verification

- A combined checkpoint passed **3,466 assertions in 92 blocks**, with zero
  failures, errors, warnings or skips. Later result/context/cleanup changes
  passed their focused suites and frozen scientific/API fixtures. Installed
  production snapshot `0067322` also passed its full package suite.
- Snapshot `0067322` passed `R CMD check --no-manual --no-vignettes` with zero
  errors/notes and the same two omitted-vignette-output warnings as the original.
  This is not a clean full-vignette check.
- Direct native component checks passed **112 assertions**; frozen native
  fixtures and exact source-oracle sampler checks remain unchanged.
- Original versus installed `f4aaab4`: all **30 actual pipeline script stages**
  passed and complete ordered payloads were exactly equal, including nuclear
  nontrial and zero-pass artifacts.
- Full standard and nuclear Nextflow workflows completed for both builds, with
  exactly equal ordered payloads and actual pod maps. Qualification uses
  Nextflow 25.10.4/Java 17/local executor, `array=null`, no containers and a fixed
  startup seed. Historical containers/default array profiles are not qualified.
- Selected unchanged book cells: **14/15 completed**, with exactly equal
  observations. Both packages fail the same CellRanger cell because two installed
  barcode files are empty. Parse and ondisc examples pass. Runners preserve the
  failed data prerequisite and return nonzero.
- The final production snapshot is `28df390`; resource measurements use its
  isolated installed library. Final qualification reports will identify the
  exact source they execute rather than relabel older runs.

Initial hosted RF0 run `38917674952` failed despite unchanged production:
Linux/Windows fitting and warning observations differ from the saved Mac
environment, and R CMD check hides LinkingTo-only BH from test-time compilation.
The new harness clean-builds verified original source in a separate process for
a paired same-environment comparison. Debug and release pairings passed locally;
frozen fixtures/API/maps/errors and exact native oracle checks are retained
without loosening tolerance. BH is also declared in Suggests for test compilation.

Hosted run `38920183528` is in progress. Its ondisc-free lane passed; Windows
fails original-source archive hash verification before candidate execution.
Archive handling and the review's environment-fingerprint P2 (compiler/config/
flags and RNG kind) are being corrected. Hosted qualification remains open.

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

## RF0 checkpoint verification

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
pipeline reports for machine-readable RF0 evidence. Current hosted results are
described separately above.

The maintainability baseline is frozen in `MAINTAINABILITY.md` and
`results/maintainability-baseline/`: fixed change-tracing exercises and a seven-
operation testability inventory, reviewed ownership sites, plus versioned R/C++
complexity/nesting and candidate-site inventories. At RF0 the unchanged candidate had
identical detailed metrics and source hashes. Human comprehension timing,
cognitive complexity, token clone coverage and a complete dynamic dependency
graph are explicitly unmeasured; there is no composite or monotonicity target.

## Independent review dispositions

The integrated production source review found no P0/P1 regression and confirmed
the real caller migrations, unchanged API/schema/samplers and removal of the
superseded factory. Its material P2 concerns stale runtime-oracle reuse after a
compiler/configuration/flags or RNG-kind change; that validation fix and its
review are pending above. Source review does not establish runtime/resource gates.

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

## Gates and limitations that remain distinct

- Hosted Linux/macOS/Windows CI and the explicit ondisc-free lane are running;
  their final results remain a gate. The backed lane must fail when ondisc is
  missing, not silently skip.
- The qualified full local Nextflow profile does not establish historical
  0.10.2 containers, nf-test, hosted executor profiles or the archived default
  dynamic-array configuration. Docker has no running daemon and nf-test is absent.
- Uncached book coverage is selected cells, not every chapter. The shared
  CellRanger example-data failure and external rd7 inputs remain explicit.
- Storage edge cases, reordered requests, discovery→calibration→discovery and
  deeper native indexing now have coverage; the legacy intercept-only NT and
  overlapping NT-local-schema limitations remain characterized, not repaired.
- Resource runs now measure compact, actual Papalexi, difficult synthetic,
  sharing-sensitive CRT and a backed assay above a declared 256 MiB dense
  materialization budget, with sampled fork/process-tree RSS. Paired results and
  threshold dispositions remain pending. The host has 48 GiB RAM; this exercise
  is not an OS-enforced out-of-memory or full-transcriptome Replogle claim.
- API/scientific/default changes, mixture repair, replacement QC, accelerated
  fits, SPA, Wilcoxon, backed CRT and new storage backends remain separate work.
