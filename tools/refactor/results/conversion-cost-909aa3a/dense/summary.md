# Dense general-conversion cost

Acceptance remains subject to the approved gate/exception review. Short timing gates are inconclusive.

Original: `21f9ea098b69c444a884a49bc254de711968e3b5`; candidate: `909aa3ad76ee21e7f9ab83a63c645bc0fb6a43c7`.
Counts: warmups=8, measured_single_conversions=32, measured_batch_children=0, profiles=8.

| Case | View | Batch | Original median s | Candidate median s | Time gate | OS peak delta MiB | OS gate | Sampled phase gate |
| --- | --- | ---: | ---: | ---: | --- | ---: | --- | --- |
| 1500000-dense | row | 1 | 0.008467 | 0.009526 | inconclusive | 2.516 | True | True |
| 1500000-dense | column | 1 | 0.008026 | 0.01008 | inconclusive | 10.109 | True | True |
| 6e+06-dense | row | 1 | 0.021435 | 0.075177 | inconclusive | 50.562 | True | True |
| 6e+06-dense | column | 1 | 0.018471 | 0.069901 | inconclusive | 49.453 | True | True |

Gate crossings: 0. No correctness-cost exception is automatically granted.

RSS is in bytes in JSON; sampled phase peaks can miss short peaks. Rprofmem records allocation traffic, not isolated live sorting workspace.
The full JSON retains every warmup, trial, profile, input mutation/output equality check, build/archive hash and threshold comparison.
Integrated workflow, backed, worker-memory/cache and consumer gates remain separate.
