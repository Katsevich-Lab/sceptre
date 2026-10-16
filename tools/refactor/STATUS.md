# Refactor status — October 16, 2026

**Active, incomplete; further cleanup implementation paused for plan review.**
Integration branch: `codex/refactor`, current safe production checkpoint
`c7d62502757f2a17741e2ef18b9e4d006bf063dc`. Nothing has been merged to main
or released. Draft [PR 234](https://github.com/Katsevich-Lab/sceptre/pull/234)
is the review surface, not a readiness declaration.

## Current checkpoint and next dependency

The initial RF1–RF8 responsibility split is integrated: canonical maps, storage
adapter, QC policy, fit/score/cache modules, separate permutation/CRT plans,
native inference components, result finalization and facade state ownership.
All six private ondisc calls have one storage owner. The public export/formal
inventory and 56-slot saved-object facade remain unchanged. PR233's new guide
cache is excluded; existing response-cache validity was corrected for actual
source design/control changes, with warm/cold regressions.

That work and passing compatibility fixtures did **not** establish whole-source
cleanliness. Eugene requires an inventory and reviewed plan before further
cleanup. [CLEANUP_PLAN.md](CLEANUP_PLAN.md) now records every handwritten module,
actual bounded callers, necessary compatibility, dead implementation, known
behavior bugs and a concrete migration queue. Coverage is 34 handwritten R
modules, 14 native implementations and three headers, with source hashes and
per-function evidence. Generated bindings are separate. Static inventory and
contracts are not proof of comprehension or absence of bugs.

Independent senior review of the first draft requested amendments. No P0;
the blocking P1 is a calibration-group bound against **post-QC active NT groups**,
not guide metadata. Exact per-entrypoint safety rules, context/factory/intercept
consistency, private-symbol/sampler dispositions, behavior ledger and measured
ownership-copy resource requirements are being incorporated before re-review.
No further cleanup implementation resumes until that plan review is recorded.

A separate source-derived P1 affects high-MOI complement calibration: summed
NT-group sizes can exceed the cell universe and reach out-of-bounds Fisher–Yates
access. The first implementation unit will cap that valid upper bound and add
precise allocation/index prerequisites without changing formerly safe draws.
The unsafe baseline case has not been executed.

Safe retained work:

- c7d6250 adds 28 core contracts with unchanged production expressions; text
  review found a source-design-order description that still needs correction.
- The canonical sparse-count/nonmutating triplet WIP is preserved in a stash
  and `codex/refactor-canonical-counts-checkpoint` at tree `8f1b8f4`. It is not
  integrated or qualified as one combined unit.
- A 39-function legacy contract draft is isolated at a2e66f0; pointer naming at
  13ff940. Neither is accepted merely because parsed code is unchanged.
- Driver/worker cleanup is a proposal only. No temporary compatibility aliases
  should be added for private names with no actual consumer.

## Completed qualification, pinned to earlier source

These are real completed gates for the named checkpoints. They do not qualify
future cleanup edits.

- Installed **90664a7**: **3,522 checks in 101 blocks**, zero failures/errors/
  warnings/skips; **112 native component checks** passed. Production expressions
  match 28df390. `R CMD check --no-manual --no-vignettes` has 0 errors/0 notes and
  the same two omitted-vignette warnings as original. This is not a full
  vignette check. [Reports](results/qualification-90664a7/qualification-summary.json).
- Hosted **c635263**, [run 38925615088](https://github.com/Katsevich-Lab/sceptre/actions/runs/38925615088):
  **all five jobs succeeded**: Windows/macOS/Linux release, Linux oldrel-1 and
  explicit in-memory/no-ondisc. The backed lane requires source ondisc 1.3.5.
  [Recorded job/step results](results/ci-c635263.json).
- Verified original-source runtime oracle uses unchanged 21f9ea0 source and
  frozen inputs, forced fresh native compilation, same resolved toolchain/
  flags/SDK/dependencies/RNG and exact full sampler banks. Environment changes
  cannot reuse an inadequately keyed oracle. Windows temporary paths are
  canonicalized for the dependency's tilde handling. No tolerance was loosened.

## Actual consumers at 28df390

Both builds are verified fresh matched-O2 installs. See
[qualification comparison](results/consumer-qualified-28df390/qualification-comparison.json)
and its retained provenance/logs/payloads.

- **30/30 actual pipeline script stages** completed for both packages with
  exactly equal ordered payloads, identities, decisions, classes and schemas.
- Full standard and nuclear Nextflow runs completed. Ordered results match under
  a shared actual-arrival/pod plan using local serialization **and fair=true**.
  Two-task parallel and serialization-only controls completed but failed exact
  gather-order equality: baseline factor/tied-row order depends on arrival.
  All five raw worker-pod payloads were exact. Failed controls remain recorded;
  outputs were never sorted to conceal the difference.
- Runtime: Nextflow 25.10.4/Java 17, arrays/containers/conda disabled and fixed R
  startup seed. Default dynamic arrays, historical containers and hosted
  scheduler profiles are unqualified.
- Selected uncached book cells: **14/15 completed exactly** for both builds.
  The same CellRanger example fails because two installed barcode files are
  empty; runners return nonzero. This is not whole-book success.

## Resources at 28df390

[Acceptance record](results/resource-release/acceptance.json): 160 primary fresh
processes, 20 settings per build across compact, real Papalexi, difficult sparse,
sharing-sensitive CRT and backed cases, cold/warm and 1/2 workers. Matching clean
-O2 builds/dependencies/configurations; all 80 paired scientific hashes exact.
All 11 meaningful elapsed gates passed (maximum median ratio 1.017).

One original Papalexi cold/two-worker aggregate-RSS outlier crossed the review
budget; main and individual-worker peaks passed. Three paired traced repeats
passed without reproducing it; the exact cause remains unresolved. The outlier
and failed build-mismatched earlier batch remain recorded. Aggregate RSS double-
counts shared fork pages and sampled peaks are lower bounds, not PSS. Host RAM
is 48 GiB; the 512 MiB dense-counterfactual backed case exceeds a declared 256 MiB
materialization budget, not an OS-enforced cap. This is not full Replogle scale.

Twelve separate allocation/fit/bank diagnostics passed. These historical cases
do not measure the pending owned-triplet conversion. That path needs meaningful
nnz/sorting-workspace/peak-RSS measurement and an explicit correctness-cost
exception review if the general gate is exceeded.

## Remaining completion criteria

1. Resolve and independently approve the cleanup addendum and exact findings/
   symbol/behavior dispositions; then implement coherent reviewed CL1–CL6 units.
2. Complete actual caller migrations, remove unused modes/helpers/fallbacks,
   correct coordinate and ownership debt and review contracts against code.
3. Requalify a final immutable build: installed suite/native checks, paired
   package checks, five hosted lanes, affected real consumers and resource gates.
4. Refresh final-source maintainability evidence and obtain independent source/
   deletion review. Existing metrics are historical, not a cleanliness score.

New methods, SPA/accelerated fits/Wilcoxon, replacement QC, new backed CRT or
aggregation, Bioconductor adapters and public-object replacement remain deferred.
Known plotting/fallback/backend limitations will have explicit dispositions;
completion cannot be described as absence of every possible defect.

The goal remains active. Routine implementation/check failures are resolved
without stopping. Important unresolved design/scope/scientific choices or
unexpected blockers are the stopping conditions. Main merge/release require a
separate explicit instruction after reviewable completion.
