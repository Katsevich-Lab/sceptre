# Refactor status — October 16, 2026

**Active, incomplete; closing known findings before one combined source freeze.**
Integration branch: `codex/refactor`. The count allocation remedy, bounded legacy
units, pair/getter cleanup, sparse vector aggregation and separate dense-entry
repair are integrated. Native retirement is source-reviewed and integrated;
the native return-lifetime repair is integrated. Final combined acceptance remains open.
Resource and combined installed-source acceptance remain open. Nothing has been merged to main
or released. Draft [PR 234](https://github.com/Katsevich-Lab/sceptre/pull/234)
is the review surface, not a readiness declaration.

## Bounded finish

The fresh independent whole-package review is retained at
[cc01683](cleanup/implementation-reviews/sceptre-whole-package-review-cc01683.md),
with per-file/flow coverage and its bounded reproduction. It endorses the R
responsibility split without a larger framework, but identifies the inherited
diagnostic-NA plotting and native return-allocation lifetime defects. The
separate [B8/B9 dispositions](cleanup/behavior-dispositions.md) bound those fixes.
The plotting/view unit passes 99 focused plotting and 88 contract checks, with
independent source review closed. The native fix passes 959 controlled cases
(12 successes, 858 raw R errors, 89 native allocation errors), a separate safe
negative control, and 112 component expectations. Scoped source/probe review
closed; the [complete native evidence](results/native-result-lifetime-6dbd1d0/README.md)
retains every attempt and the recovered transformed source. Global/pre-leaf
Rcpp precious-object cleanup is outside this declared native-buffer guarantee.

Finish these known findings and the already prepared contract/view corrections,
then select one explicit immutable production checkpoint. Further nonblocking
improvements go to follow-up work; reopen production only for a required gate
failure or material finding. No percentage estimate describes the remaining time.

The final acceptance record must bind the same checkpoint to all rows below.
Earlier passing counts remain historical evidence, not a combined qualification.

| Required final gate | Evidence needed |
| --- | --- |
| Immutable source/build | Archive, source blobs, fresh native compilation, DLL and matched original toolchain/dependencies |
| Installed package and native components | Full suite with required backed coverage, API/saved/workflow/index/cache/RNG contracts, independent original oracle and full active sampler banks |
| Native error ownership | Count, draw, dispersion and new result-return controlled probes; explicit scope and negative controls |
| Actual consumers | All 30 script stages; standard/nuclear Nextflow under the disclosed shared arrival plan; selected book cells with the known baseline failure retained |
| Resources | Fresh complete-workflow comparison and call/allocation diagnostics; source-invariant carry-forward of separately qualified B6 and accepted triplet/dense paths, with material deltas reassessed |
| Maintainability/review | Current measurements with fixed denominators, four repeated source traces, complete architecture/native-contract review and fresh independent follow-up |
| Bioconductor checks | Actual paired BiocCheck on original/frozen candidate with matched supported R/Bioc/checker versions; full logs and inherited/new/submission-only dispositions |
| Hosted checks | All five jobs successful on the final code-bearing checkpoint; retained failed runs remain visible |

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

- c7d6250 adds 28 core contracts with unchanged production expressions; CL3
  corrected the source-design-order description found inaccurate in text review.
- The original canonical-count WIP remains recoverable at tree `8f1b8f4`; its
  reworked coherent implementation is now integrated as `c84e85f` below.
- A 39-function legacy contract draft is isolated at a2e66f0; pointer naming at
  13ff940. Neither is accepted merely because parsed code is unchanged.
- Driver/worker implementation `627e7e4` is integrated as `0b6ae38` after
  focused immutable qualification and independent source review. No temporary
  compatibility aliases were added for private names with no actual consumer.

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

CL2 (`ab477759`, integrated as `c84e85f`) passed 1,059 focused expectations.
Its first source review closed, but a subsequent allocation-path audit found a
material P2: raw R allocation errors can bypass newly owned C++ coordinate and
scratch destructors. The explicit review amendment supersedes that part of the
first closure. The narrow phase-order remedy at `8841efa` (integrated as
`20945f0`) has independent source closure, 1,065 installed expectations in 50
blocks and 214 controlled error/success cases plus an unsafe-order negative
control. The first complete exact-source triplet study is retained at 8841efa.
All eight returned outputs match exactly and candidate inputs remain unchanged;
all four large five-conversion timing gates pass. All four large settings exceed
both memory gates (whole-process increases of 111–125 MiB), while small settings
pass. The [unchanged complete experiment](results/triplet-experiment-8841efa/summary.md)
is not resource acceptance. A source-static
copy reduction at da6a079 passed 1,095 installed checks. The revised 909aa3a
study completed: exactness and all meaningful timing gates pass; all small and
logical memory gates pass. Two eight-million-entry numeric cases retain
121.938/123.141 MiB whole-child increases and cross the generic gate. Eugene
explicitly accepted these measured bounded correctness costs in
[decision D2026-10-16-02](cleanup/decisions.md). The crossings remain failures
under the original gate, with an accepted exception; no optimization/rerun is
needed solely to make these unchanged cases green. Evidence and the separate
dense-cost study are [retained in full](results/conversion-cost-909aa3a/README.md).
Mechanism review and final integrated source/resource qualification still apply.
The indexing-test handoff to CL3 is complete. CL3 (`627e7e4`) passed 1,006
expectations in 52 tests against a fresh same-environment original oracle;
independent source review closed and integration is complete. This does not
qualify the combined final source or downstream consumers/resources.

B1 (`226a5b4`) has 46 candidate-source import regressions and independent source
review. Its owned in-memory Cell Ranger/Parse constructors preserve sorted
coordinates/zeros/IDs, map actual modality rows and reject repeated coordinates
per Eugene's choice. Installed-source qualification is pending; sparse sorting
workspace is unmeasured. The externally delegated ondisc file parser is unchanged
and unqualified for this policy. See the explicit boundary in decisions.md.
B6's initial Mac-only characterization and 80-child resource experiment at
88ed37c remain preserved in [the original evidence](results/vector-aggregation-88ed37c/README.md).
Hosted Linux CI exposed a real compatibility defect: dense column sums use an
extended accumulator there while the replacement sparse sum uses doubles.
Those initial results therefore do not establish cross-platform acceptance.
The repair keeps sparse row preparation but sums complete selected rows in
bounded dense column blocks, preserving the current platform's original
accumulator. Tests now execute the frozen original on the same platform,
including all 28 immutable fixtures and 200 seeded duplicate/rounding cases.
The optimized 783a250 build passes 613 installed checks and 228 exact original
comparisons. Its complete 80-child resource study passes four exact-output/RSS
pairs and all three meaningful timing gates; actual-data timing remains
inconclusive below one second. Independent source/resource review closed.
[All passing and failed attempts](results/aggregation-qualified-783a250/README.md)
and [the superseding amendment](cleanup/aggregation-precision-amendment.md)
remain preserved. Hosted 3d9fafc numerical tests pass in all five lanes; standard
checks fail only the new tests' undeclared digest dependency, corrected in
Suggests at e933db5. Final combined installed/hosted qualification remains open.
Dense output/transpose workspace, metadata multiplicity and backed capability
are unchanged.
Dispersion taxonomy has 107 exact
original-state comparisons. The two-type narrowing at `ef0c77b` passed 107
actual-source numerical comparisons and 14 paired controlled MLE/MM failure
cases, including partial-state and unrelated failure propagation, plus
independent source/proof review. All five hosted native policy probes passed at
c47e7ae; final combined installed-source evidence remains required.
File-lifetime `abf7245` has 12 source-evaluated
regressions and independent review; installed-source qualification is pending.

B2/O1 at 2f40f4f has 61 independent candidate-source checks and six exact complete
old/new trans outputs. It repairs unsorted cis selection and wide midpoint
arithmetic, and owns trans exclusion tables. Getter 09a704c has 24 independent
source checks and six exact original comparisons, replacing the final vector
offset caller with owned R arithmetic. Independent source review found no
material issue; final combined installed-source qualification remains open. Their initial test-only
fixture/attribute mismatches are retained in the unit evidence.

The legacy plotting source review and its held tiny old/new runtime comparison
are complete: theme, layer data and appearance settings match in both MOI
fixtures. A separate newly found dense symmetric/near-symmetric conversion loss
has a reproduced P2 disposition: use supported general conversion before
sparsification to retain all supplied values, then measure its dense-input-only
copy cost. Source bc53bc1 (integrated as 78636bc) has independent review and
1,125 installed checks in 53 blocks, including the separate dense regressions.
Its separate bounded dense-cost study passes all four OS memory gates; timing
is explicitly inconclusive (8–75ms), and missing sampled-phase observations are
counted rather than converted into passes. Independent cost review closed at
909aa3a. It is not folded into the numeric-triplet exception.

Hosted c47e7ae [passed all five jobs](results/ci-c47e7ae.json), including
draw-bank ownership, count ownership and dispersion-policy probes in every
lane. Windows progress-byte expectations preserve CRLF; all four standard jobs
use the verified corrected ondisc source pin. This is a prior checkpoint,
not qualification of the newer combined source. The retained-R source/contract
audit closed at 750b1b5, with exact mixture prose closure at 88a613a; native
retirement, full native documentation and final architecture review remain open.

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
5. Review every retained handwritten C++ implementation/header for useful
   maintainer documentation, including numerical meaning, cell/index/draw views,
   reuse, ownership/allocation ordering and fallback state. Audit the whole
   ARCHITECTURE.md against final names, actual callers and compatibility
   boundaries. Eugene explicitly requires both in source/documentation coverage;
   passing tests or metric improvements cannot substitute for this review.

New methods, SPA/accelerated fits/Wilcoxon, replacement QC, new backed CRT or
aggregation, Bioconductor adapters and public-object replacement remain deferred.
Known plotting/fallback/backend limitations will have explicit dispositions;
completion cannot be described as absence of every possible defect.

The goal remains active. Routine implementation/check failures are resolved
without stopping. Important unresolved design/scope/scientific choices or
unexpected blockers are the stopping conditions. Main merge/release require a
separate explicit instruction after reviewable completion.
