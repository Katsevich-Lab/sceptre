# Refactor status — October 16, 2026

**Active, incomplete; safety qualified, count/context cleanup under way.**
Integration branch: `codex/refactor`, current safe production checkpoint
`489ca7c` (production identical to qualified safety commit `5ac0fca`). Nothing has been merged to main
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

Independent senior review requested and then verified the plan amendments.
The [signed amendment review](cleanup/senior-review-amendment.md) approves the
bounded implementation with no remaining P0/P1/material P2 planning gaps. The
complete reviewed artifact bytes are preserved at 08e0260. Exact post-QC native
bounds, context/factory/intercept consistency, symbol/sampler dispositions,
behavior ledger and measured ownership-cost requirements are now settled.
CL1 preceded broad restructuring and is now qualified locally below. The plan
review did not itself qualify code or approve a resource exception.

The audit also identified a high-MOI complement-calibration P1: summed NT-group
sizes could exceed the cell universe and reach out-of-bounds Fisher–Yates access.
CL1 caps that upper bound and checks precise allocation/index prerequisites
without changing formerly safe draws. The unsafe baseline was not executed.

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

## CL1 safety integrated

Reviewed safety commit 5ac0fca is integrated as 489ca7c. Its exact immutable
release build passed 3,631 installed checks in 110 blocks, 344 focused checks,
112 native component checks and 12 controlled R-boundary ownership cases, with
zero failures/errors/warnings/skips. The fresh same-environment original oracle
preserves all frozen workflow/API/cache/saved/serial-fork/backed cases. See
[qualification and provenance](results/qualification-cl1-5ac0fca/qualification-summary.json).
The initial late-registration P2 was corrected: a protected finalizing NULL
external pointer precedes native allocation; only nonallocating ownership
transfer follows. Independent source review closed the finding. Failed stale-
oracle and reporter-setup attempts remain identified separately.

CL2 native/count ownership and the coupled CL3 context/execution migration now
have separate worktree/file owners. The indexing test has an explicit sequential
handoff; no shared production file has concurrent writers. B1 repeated-file
coordinates will be rejected per Eugene's recorded choice. B6 characterization
has 28 exact helper and 200 exact rounding matches, plus public import probes;
resource qualification remains required. Dispersion taxonomy has 107 exact
original-state comparisons and a source-reviewed two-type narrowing path;
controlled failure propagation and platform facts remain required.

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

1. Finish the independently approved CL2–CL6 units after accepted CL1 safety;
   keep conditional characterization/resource gates and named deferrals explicit.
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
