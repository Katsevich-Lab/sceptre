# B6 sparse vector aggregation independent review

Exact candidate `88ed37cd199c642edd0ba0c008ab507bab01d10f`, base `7df5ad1b164b00ae45f91a0ded3085c45764a8ca`, original numerical/resource reference `21f9ea098b69c444a884a49bc254de711968e3b5`. Reviewed both ordered commits `3a013c6` then `88ed37c`, complete changed source/tests/tooling, declared resource protocol and all completed recorded trials. Read-only; no R execution, build, test, benchmark or source edit.

**Outcome: no P0/P1/material P2 found. The bounded B6 unit meets its exact-output, primary/meaningful timing and all four RSS gates on the recorded environment.** The real-data timing remains inconclusive, and combined integration qualification is still required. This is not a claim that all numerical inputs, Matrix versions or package workloads were exhaustively proved equivalent.

## Source and correctness

The implementation makes one focused representation change: accepted sparse input becomes row-accessible once through public Matrix coercion instead of realizing the whole assay densely. It retains the existing first-occurrence vector IDs, dplyr row selector (including repeated metadata selections), ordered selected guide rows, Matrix::colSums, numeric vapply assembly and transpose. Non-sparse coercion/dispatch remains the old as.matrix route. The required vector-by-cell output and assembly/transpose workspace remain dense, explicitly documented. No replacement summation algorithm, compensated rounding, count floor or custom sparse duplicate policy was introduced.

The source deliberately avoids storage's stored-position conversion, whose logical FALSE-to-one rule would be wrong for count aggregation. Matrix handles accepted numeric/logical sparse duplicates and FALSE values; tests independently distinguish coordinate duplicate aggregation from repeated metadata contributions. Caller inputs remain borrowed/read-only through Matrix coercion/subsetting. There is no new cache, RNG, backing-file route, per-cell group index bank or full sparse-input dense realization.

The independent fixtures cover numeric/logical dense and C/R/T inputs, valid duplicate triplets, explicit zero/FALSE, repeated metadata, first-seen vector/source-cell order, dense double class/dimnames, public import metadata and serialized input nonmutation. The rounding fixture specifically preserves the observed original double result at 2^53 and for .1+.2+.3. Installed characterization uses separate selected namespaces for 28 frozen plus 200 generated duplicate/rounding inputs; its runner checks exact output equality and equal generated input digests. Historical fixtures are hash-pinned and not regenerated. These records support compatibility on the tested cases; they are not substituted for source reasoning.

## Resource protocol integrity

The protocol was committed before the run. Actual execution contains **80 fresh children**: 8 correctness, 32 single-call timing, 32 20-call timing-batch and 8 separate allocation-profile children. Each mode has one warm-up pair plus three measured pairs where specified; build order alternates as declared. All original records, phase samples, time-l output, profiles and attempt hashes are retained. The first sandbox ps failure is honestly preserved as a failed smoke, followed by a successful narrow escalated smoke; it is not relabeled as a resource trial. The sole full run has no failure/safety stop or hidden threshold repeat.

Input generation is sparse, with a 3,200-by-24,000 primary assay (0.5% stored density), repeated numeric triplets, stored-FALSE logical input and actual Papalexi guide counts under explicitly synthetic four-guide grouping. Exact complete result RDS byte hashes were independently recomputed and match for all four reference/candidate pairs. I checked all 80 recorded environment fingerprints: one common dependency/numerical/thread environment, no warnings, and all child return codes/safety flags successful. Minimum observed available-memory proxy was 19.24 GiB, above the declared 8 GiB launch floor.

Timers surround only the helper; input I/O, explicit pre-call/inter-call GC, validation and monitor acknowledgements sit outside. Automatic GC within the helper is measured. Batch children independently reload inputs, keep one input/output live at a time and sum the 20 conversion intervals; batch RSS is not substituted for the primary single-call memory measure. Profiles are separate and excluded from timing/RSS gates.

OS whole-child maximum RSS is authoritative; 20ms actual-R-PID sampling provides acknowledged pre/post boundaries and active observations, with its brief-peak limitation disclosed. All measured single-call trials had active samples. R heap limits and conservative sampled 2 GiB stop are explicitly not hard native-memory guarantees. No deliberate OOM or concurrent-heavy-work evidence is claimed.

| Case | Meaningful batch median original -> candidate | Candidate/original | Maximum measured single-call OS RSS original -> candidate | Gate disposition |
| --- | --- | --- | --- | --- |
| large_numeric_C | 20.912s -> 12.353s | 0.591 | 1374.1 -> 613.6 MiB | pass; RSS pass |
| medium_numeric_T_repeats | 2.043s -> 1.686s | 0.825 | 326.5 -> 279.4 MiB | pass; RSS pass |
| medium_logical_C | 1.982s -> 1.719s | 0.867 | 296.1 -> 280.5 MiB | pass; RSS pass |
| actual_papalexi_guides | 0.673s -> 0.467s | 0.694 | 293.9 -> 271.7 MiB | inconclusive (<1s); RSS pass |

All three synthetic batch medians exceed the declared one-second meaningful threshold in both builds; all satisfy the <=10% time increase gate. Papalexi's 0.673/0.467s batch remains explicitly inconclusive, so there is no real-data speed claim. All four OS and observed active RSS gates pass. The primary allocation profile identifies the original 614,400,048-byte full-input allocation; candidate has none at that size and its largest R allocation is 153,600,048 bytes, the retained dense output/assembly scale. That distinguishes removal of input densification from a fictitious sparse-output claim. Traces omit small/native allocations as disclosed; RSS still includes process memory.

I independently checked both source archives' recorded production Git blobs, current installed DLL hashes, identical release compiler prefixes and the recorded harness/worker SHA256 against the exact candidate commit. These checks and complete output hashes are in `/private/tmp/sceptre-b6-independent-evidence-check-88ed37c.json`. No special toolchain/faster flags or candidate-generated original oracle was substituted. Durable compressed all-trial and allocation evidence is self-contained; large input/output files additionally remain pinned in the original attempt directory.

## Remaining qualification boundary

Integrate both ordered source/tool commits as one reviewed unit, preserving root's later independent comment changes. Run final combined installed/public/API/saved-workflow gates on the actual integration source; this candidate branch predates later count ownership, dense repair and native retirement work. The supported external ondisc aggregation route is unchanged. The optional actual-data timing cannot become a pass by quoting its attractive subsecond ratio. No further B6 redesign or speculative feature is requested by this review.

## Exact source and evidence hashes

- `R/import_functs.R`: `b875b45725e5b39e5099417980d9ca6c99db5fd82d750c17b20a25eed3f8b0d4`
- `tests/testthat/test-refactor-vector-aggregation.R`: `cb1058338888d5e2f028823adc383894c0b7f73632ed1bb96c262116d99c330e`
- `tools/refactor/aggregation-characterization.R`: `129a9287e7b63ee393a6d9ed52ef8c08f6affe3741eb81894b86a7ee79e340d2`
- `tools/refactor/aggregation-data.R`: `5e71bacd78f19151798da0a3f5bfd33275ad485ba0f86796c71647ce83bad3c6`
- `tools/refactor/aggregation-worker.R`: `ae84c40b7f056fa09b108273c5747b30625f5b5337864d1a118b1dc4c0310c10`
- `tools/refactor/aggregation-benchmark.py`: `e3c8e4cad9a9a49bf939eaeb3d753c4bffeccb8562189da71f2ec8bcac1059cf`
- `tools/refactor/cleanup/vector-aggregation-resource-plan.md`: `80b1a297854cb9e9549b7375e6e5f8ba1cf474e82e34be2d78abcf81cc389b67`
- `/private/tmp/sceptre-b6-resource-88ed37c-run1/summary.json`: `711f099d0d5e128b626768eefa48fea712a2bf858322c14ef0529a8326e40417`
- `/private/tmp/sceptre-b6-evidence-88ed37c/report.json`: `3eb2a58ecc131aac97a191a551ac8deeb6bb318352f354ec7aa17fe459e3203a`
- `/private/tmp/sceptre-b6-evidence-88ed37c/all-trials.json.gz`: `c315d16cd16d2f32a9a3782d2692f85cf58d7931e99d8591a4360773ca3580a8`
- `/private/tmp/sceptre-b6-evidence-88ed37c/allocation-traces.json.gz`: `3dcc55c33f3c0d09e15c227f337dcd436ac9469d043ae474d8a4b68a7b8d5f24`
- `/private/tmp/sceptre-b6-evidence-88ed37c/attempt-file-manifest.json`: `3bd4c21a8d3289accd42b3ec2b0f10386d4cbf0f920f93a02060ec5ac591bf45`
- `/private/tmp/sceptre-b6-independent-evidence-check-88ed37c.json`: `7565dc352f755e2570c64f8c913f40977d9e7ad59c3846a79a991742e3154e0b`

## Durable evidence checkpoint readback

Final evidence tip `f070b1ca63c09b30d31347ae2b0d3d60fbf6f915` has no difference in R/src/DESCRIPTION/NAMESPACE from measured `88ed37c`. Its committed README was reviewed and correctly preserves the actual-data timing limitation, dense result/transpose costs, initial failed smoke and final integration boundary. The four principal report/raw-evidence files exactly match the independently inspected temporary evidence. Every committed artifact-manifest size/SHA256 was independently checked. The 80 full trial records and eight allocation traces are retained in compressed JSON, not replaced by successful summaries.
