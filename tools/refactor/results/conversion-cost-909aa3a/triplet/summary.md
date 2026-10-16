# Triplet conversion resource measurements

Acceptance remains subject to the approved gate/exception review. Short timing gates are inconclusive.

Original: `21f9ea098b69c444a884a49bc254de711968e3b5`; candidate: `909aa3ad76ee21e7f9ab83a63c645bc0fb6a43c7`.
Counts: warmups=16, measured_single_conversions=64, measured_batch_children=32, profiles=16.

| Case | View | Batch | Original median s | Candidate median s | Time gate | OS peak delta MiB | OS gate | Sampled phase gate |
| --- | --- | ---: | ---: | ---: | --- | ---: | --- | --- |
| 1000000-numeric | row | 1 | 0.036373 | 0.038129 | inconclusive | 15.438 | True | True |
| 1000000-numeric | column | 1 | 0.048831 | 0.046355 | inconclusive | 15.078 | True | True |
| 1000000-logical | row | 1 | 0.038121 | 0.039487 | inconclusive | 2.422 | True | True |
| 1000000-logical | column | 1 | 0.039249 | 0.039285 | inconclusive | 4.266 | True | True |
| 8000000-numeric | row | 1 | 0.48474 | 0.537422 | inconclusive | 121.938 | False | False |
| 8000000-numeric | row | 5 | 2.212975 | 2.396044 | True | n/a | batch excluded | batch excluded |
| 8000000-numeric | column | 1 | 0.567248 | 0.633317 | inconclusive | 123.141 | False | False |
| 8000000-numeric | column | 5 | 2.64785 | 2.831311 | True | n/a | batch excluded | batch excluded |
| 8000000-logical | row | 1 | 0.546809 | 0.599274 | inconclusive | 61.891 | True | True |
| 8000000-logical | row | 5 | 2.304938 | 2.478995 | True | n/a | batch excluded | batch excluded |
| 8000000-logical | column | 1 | 0.631485 | 0.67928 | inconclusive | 58.188 | True | True |
| 8000000-logical | column | 5 | 2.828742 | 2.952814 | True | n/a | batch excluded | batch excluded |

Gate crossings: 2. No correctness-cost exception is automatically granted.

RSS is in bytes in JSON; sampled phase peaks can miss short peaks. Rprofmem records allocation traffic, not isolated live sorting workspace.
The full JSON retains every warmup, trial, profile, input mutation/output equality check, build/archive hash and threshold comparison.
Integrated workflow, backed, worker-memory/cache and consumer gates remain separate.
