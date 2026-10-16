# Bounded dense-conversion cost independent review

Reviewed exact tooling/recipe commit `909aa3ad76ee21e7f9ab83a63c645bc0fb6a43c7`, its diff, frozen protocol/tool copies and `/private/tmp/sceptre-dense-cost-909aa3a/summary.json` plus small recorded child/monitor metadata. Read-only during the triplet quiet window. No R process, build, test or benchmark; no large RDS/assay read, rehash or deserialization.

**Outcome: no P0/P1/material P2 in this bounded source/tooling/evidence review.** Four exact/RSS comparisons complete successfully. All timing comparisons remain inconclusive as predeclared; no timing pass or automatic correctness-cost exception is granted. Final integrated qualification and the independent triplet/resource gates remain separate.

## What was measured

The actual reference is immutable original `21f9ea098b69c444a884a49bc254de711968e3b5`, against repaired candidate `909aa3a`; this is **net dense-path cost against original**, not an isolated da6-parent-versus-repair experiment. Source review confirms the original rectangular dense branch already coerces to fresh T slots and sorts them in place through an owned list/setDT. The candidate also sorts fresh T slots, after public general-before-T coercion, with small static ownership/list-setup changes. Do not attribute the measured net time/RSS delta wholly to one line or call it an extra defensive copy of sparse columns.

The supported public general coercion adds a dense payload allocation, as previously established from Matrix source. The observed profile deltas of about 8 bytes per dense input element are consistent with that extra payload, plus small other changes, but are allocation traffic rather than a pure component or live-workspace measurement. The two rectangular fixtures intentionally avoid the baseline's symmetry data-loss bug; separate small symmetric/near-symmetric regressions, not these equal-output performance fixtures, establish the correctness repair.

Production R/src/DESCRIPTION/NAMESPACE in `909aa3a` is byte-identical to previously source-reviewed `bc53bc1`; 909 changes only resource tools/protocol. It is not the final main integration source, and its measured net comparison should retain these exact pins.

## Recipe/tooling and provenance

The small extension reuses the existing clock, worker, phase sampler, exact sparse-output comparator, build/archive/package pins and safety monitor. `--recipe dense` chooses a distinct frozen protocol/generator and requires exactly four fixture/view settings. The worker takes the declared dense nonzero count without accessing an absent dense @x slot. Triplet remains the default recipe, retains its expected settings, and now also refuses to summarize a short large-case timing series without its required batch records. No arbitrary recipe framework or production hooks were introduced.

The deterministic dense fixtures have 1.5M/6M elements (12/48MB double payloads), exactly 10% nonzero linear positions and both requested orientations. Source generator arithmetic/names match the declaration; there are no random draws or refreshed scientific goldens. Exact comparison reads class, dimensions, dimnames and every p/index/x slot; recorded caller serialization hashes are unchanged.

I checked the five tool hashes and protocol hash against exact `909aa3a` and frozen harness copies; both recorded source trees still match their recorded production Git blobs and release compile prefixes are identical. All **56 helper-conversion records** have one common R/dependency/DLL/thread environment after excluding the deliberately distinct sceptre library. Installed-package/archive/input/output hashes are present in the preserved records; large RDS hashes were not independently recomputed here, per the quiet-window scope. All four recorded complete output hashes agree and the actual comparator records full exact slots rather than only a checksum/total.

Counts are internally consistent: 8 correctness conversions plus 4 comparator children; 8 warm-ups; 32 primary measured conversions (four fresh pairs per fixture/view); 8 separate allocation profiles. That is 56 helper-conversion children and 60 conversion/comparison monitor records, with all return/safety fields successful. The fixture preparation is additional, not a silently counted conversion. Both views retain original-first/candidate-first/candidate-first/original-first measured ordering. No timing-only batch was planned or added.

Headroom checks at preparation/correctness/measurement phase start recorded approximately 19.26/19.79/19.27 GiB, above 8 GiB. The same 1536M R heap, one child at a time, one thread and 2GiB sampled/OS retrospective stop apply. This is a reclaimable-memory proxy and monitored bound, not an all-native-allocation guarantee or a claim about unrelated host activity. No failed/safety-stopped child is omitted from the successful summary.

## Observed cost and limits

| Dense elements | View | Original -> candidate median | Maximum OS RSS increase | Profile traffic increase per dense element | Disposition |
| --- | --- | --- | --- | --- | --- |
| 1,500,000 | row | 8.467 -> 9.526 ms | 2.516 MiB | 7.9837 B | RSS pass; time inconclusive |
| 1,500,000 | column | 8.026 -> 10.080 ms | 10.109 MiB | 7.9837 B | RSS pass; time inconclusive |
| 6,000,000 | row | 21.435 -> 75.177 ms | 50.562 MiB | 7.9959 B | RSS pass; time inconclusive |
| 6,000,000 | column | 18.471 -> 69.901 ms | 49.453 MiB | 7.9959 B | RSS pass; time inconclusive |

All four OS RSS increases are below the declared max(10%,64MiB) allowance; the largest is 50.562MiB. The absolute increases remain visible instead of being called zero cost. All original and candidate medians are far below the predeclared one-second threshold. In particular the attractive or adverse relative ratios are not meaningful-time passes/failures under this bounded recipe. `accepted_correctness_cost_exception` remains null.

There were **19 of 40 timing/warm-up children with no in-conversion RSS sample**, expected when an 8–75ms interval meets a 20ms sampler. The per-setting sampled-phase gates therefore summarize observed subset maxima, not complete per-trial peak coverage. All-child OS whole-process maxima remain the primary complete memory authority. Sixty-millisecond untimed pre/post padding helps observe live-input/output boundaries; it cannot recover a missed transient conversion peak. This limitation must remain explicit if sampled-phase passes are quoted.

The 12/48MB candidate profile maxima and approximately 8B-per-input-element additional traffic support the cost explanation. They do not prove that dense copying dominates all future shapes, bound the cost for arbitrarily larger dense inputs, or isolate sorting workspace. Sparse/backed paths are not densified by the repair; their own paired resources are not waived.

The recipe and summary are candid about short timings, allocation-versus-RSS, source pins and final integration. The durable explanatory report should call these deltas net dense-path costs versus original and preserve the missing-phase-sample limitation. No new method, private Matrix constructor, output class policy, fallback or performance waiver is requested.

## Evidence binding

- `/private/tmp/sceptre-dense-cost-909aa3a/summary.json` SHA256 `317013eb9518dcbf80e38edd00c5ba5c67f2a5b5a24e3bf598e18e9eb00116ca`
- `/private/tmp/sceptre-dense-cost-909aa3a/summary.md` SHA256 `705043116e084ccf19f1092851421001dc0e289aea55238e3bbbdaf953db6bf7`
- `/private/tmp/sceptre-dense-cost-909aa3a/plan.json` SHA256 `74da190af915285a69bb53dacde59be505170cb028e576fbdd3eaa7f9af72002`
- `909aa3ad76ee21e7f9ab83a63c645bc0fb6a43c7:tools/refactor/DENSE_CONVERSION_COST.md` SHA256 `f01eadcccb64fa87cf864bf0136263a0a379e2ed41fc2bd1cc41c4eacb4d9b94`
- `909aa3ad76ee21e7f9ab83a63c645bc0fb6a43c7:tools/refactor/dense-conversion-data.R` SHA256 `55663abd2c036c8fad3cf76e89d701e74b5e418fd1928f529347a8c9f884eb11`
- `909aa3ad76ee21e7f9ab83a63c645bc0fb6a43c7:tools/refactor/triplet-resource.py` SHA256 `286f062eef0bc52c90f0dc79d8dba35913c85066c269aa2fc28e334fa527c5c3`
- `909aa3ad76ee21e7f9ab83a63c645bc0fb6a43c7:tools/refactor/triplet-resource-worker.R` SHA256 `3ffea05e9fff4598e19e2eec69559b3074cdf57415187f797d8de57d821074bd`
