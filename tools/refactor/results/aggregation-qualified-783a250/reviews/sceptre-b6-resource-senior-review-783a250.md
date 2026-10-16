# Independent B6 optimized resource evidence review

Production pin: 783a2501295acc0b3fddab7bd2a615e545e5ae37; original pin: 21f9ea098b69c444a884a49bc254de711968e3b5. Review is read-only and source/evidence only. No R process, build, benchmark, large RDS deserialization or result re-generation was run. This is a new report: earlier 963 source review and failed resource study remain intact.

## Outcome

No P0/P1/material P2 found in this completed, predeclared four-case resource evidence. The three meaningful timing gates and all four exact-output/OS-RSS gates pass at the measured pin. Actual Papalexi timing remains explicitly inconclusive, not a pass. This closes this bounded Mac resource study; hosted cross-platform precision checks, exact final-source carry, combined installed tests and final package qualification remain separate.

## Independent checks

Read summary/plan, protocol and worker/harness source, all 80 child result/memory/time records, eight raw allocation traces, both build manifests and immutable production source snapshots. Recomputed all 51 candidate and 39 original production Git blob hashes; both installed DLL SHA256s still match. The summary plan equals saved plan; harness/worker hashes match the recorded plan. Both recorded clean -O2 builds use the same Apple clang21/SDK and dependency environment. All80 children report the same R4.6.0 arm64 Mac/Matrix1.7.5/dplyr1.2.1/data.table1.18.4/Rcpp1.1.1.1.1/BH1.90.0.1/ondisc1.3.5 environment, one data.table/BLAS/OpenMP thread and R_MAX_VSIZE=1536M.

Counts reconcile: 32 single children +32 batch20 children (each mode16 warmup/measured pairs across4 cases) +8 correctness +8 profiles =80. Each mode/case/build has one warmup and three measured fresh children, with alternating build order. All return0, no safety stops/warnings, input nonmutation true, observed OS maxima. Minimum headroom proxy19.218GiB; largest OS peak across every child1720.984MiB, below conservative2GiB stop. Every64 timing/batch child has an active RSS sample; sampling is still not an exact peak detector.

Recomputed medians and OS maxima from trials and allocation counts/bytes/maxima directly from raw Rprofmem traces. OS whole-child RSS remains authority; profiles are separate from timing/RSS gates and exclude allocations<64KiB/native malloc. Timers exclude explicit input load/GC/validation/monitor waits, include automatic GC. Batch times sum20 independent conversion intervals with one input/output live per interval, not process elapsed.

| Case | batch original/candidate seconds | ratio | single maximum OS RSS original/candidate MiB | result |
|---|---:|---:|---:|---|
| large numeric CSC |21.796/15.284|0.70123|1354.469/622.375|exact, time/RSS pass|
| numeric triplet +repeats |2.090/2.011|0.96220|326.891/300.172|exact, time/RSS pass|
| logical CSC +stored FALSE |1.952/1.917|0.98207|296.219/286.594|exact, time/RSS pass|
| actual Papalexi counts/synthetic groups |0.665/0.540|0.81203|296.000/288.062|exact/RSS pass; time inconclusive<1s|

The correctness harness compares full serialized result bytes before any timed children, not only sums. Those four comparisons are recorded true; output class/dimnames/order/type remain checked and inputs unchanged. I did not rehash the large correctness RDS outputs. The recorded separate installed log shows613 expectations/19 blocks without failures/errors/warnings/skips. The installed replay runner enforces the frozen28 input hash plus200 seeded cases, distinct installed namespace path/DLL and input nonmutation; its small candidate/reference artifacts are pinned below. The parent reports228 exact outputs; I did not re-deserialize those RDS files in this source/evidence review (whole records intentionally differ in package/DLL metadata).

## Mechanism and limits

Profiles show the original full-input dense allocation in all four cases; candidate has none at or above full-input size. Largest candidate traced allocations are the required dense result (large153600048bytes, medium9600048, actual4643344); no dense whole-assay allocation is claimed hidden by profile traffic. Each predeclared group has up to four guides, so these resource cases primarily qualify the new direct-original dense colSums fast path on small selected blocks. The >2^20-elements column-chunk path has separate exact block-boundary source/tests; its broader resource performance is not established by these four grouped workloads. Source keeps full selected rows per sum and platform dense accumulation, no partial-row accumulator or sparse colSums rounding substitute.

These are one machine/environment and four declared inputs, not a universal speed/memory claim, precision proof for every platform, biological validation, or full-package qualification. Actual guide counts have synthetic grouping. The prior963 medium-triplet timing crossing (ratio1.150) stays retained; it is not relabeled as a783 pass. Earlier88ed Mac-only sparse accumulator acceptance was superseded by the Linux precision finding. Final qualification must retain both those histories and prove the final helper carries this reviewed arithmetic/source behavior.

## Byte pins

- `/private/tmp/sceptre-b6-fastpath-783a250-resource-run1/summary.json`: `24de5f3f87703bcdb6cbb7db7d8f0b242f9ae97914b5ba0d53840b6df4103e3c`
- `/private/tmp/sceptre-b6-fastpath-783a250-resource-run1/plan.json`: `54f21ce35810187bbe7fcf75fd052ab811590c0b7a4a376af7608174fec50596`
- `/private/tmp/sceptre-b6-resource-independent-checks-783a250.json`: `48f609a973195ba6035a13de5bc78b039bf6fd6b7dd155edf7ba211b881ab305`
- `/private/tmp/sceptre-b6-fastpath-783a250/candidate/build.json`: `a8e0bf3629671cecdd3b02488b37e9a39342ca177b4576f60ea8b3440cf7652b`
- `/private/tmp/sceptre-b6-precision-9631ef2-attempt2/reference/build.json`: `1e77ac063e0426120de2fda28b3dfa3a5efb26108e350f6f871a9e502c8d93b5`
- `/private/tmp/sceptre-b6-fastpath-783a250/installed-tests.log`: `5ab157c9314077db7c5d59409d1b91424c94b53817c3cfd39f132800cb082f20`
- `/private/tmp/sceptre-b6-fastpath-783a250/replay-candidate.rds`: `488de0947d76fc850f193384f041e8b0d118d336da8d29c3ef9359a99d7451b6`
- `/private/tmp/sceptre-b6-precision-9631ef2-attempt2/replay-reference.rds`: `899224fef9fff4404534a8064c8e378a9572d95c3008d32ddf861ad4a788ea8b`
- `/private/tmp/sceptre-b6-precision-9631ef2-resource-run1/summary.json`: `b34ea54d7f1bd74edbc8e5eb2fb211bdd60a23079c71ca457c988b005f6e0a28`
