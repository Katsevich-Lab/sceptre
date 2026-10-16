# Triplet conversion resource measurements

Acceptance remains subject to the approved gate/exception review. Short timing gates are inconclusive.

Original: `21f9ea098b69c444a884a49bc254de711968e3b5`; candidate: `8841efa0a953fc545bf3b09197ab0e74e1968e1b`.
Counts: warmups=16, measured_single_conversions=64, measured_batch_children=32, profiles=16.

| Case | View | Batch | Original median s | Candidate median s | Time gate | OS peak delta MiB | OS gate | Sampled phase gate |
| --- | --- | ---: | ---: | ---: | --- | ---: | --- | --- |
| 1000000-numeric | row | 1 | 0.036226 | 0.038028 | inconclusive | 15.391 | True | True |
| 1000000-numeric | column | 1 | 0.03697 | 0.03796 | inconclusive | 19.547 | True | True |
| 1000000-logical | row | 1 | 0.038002 | 0.039719 | inconclusive | 13.844 | True | True |
| 1000000-logical | column | 1 | 0.037398 | 0.039098 | inconclusive | 16.625 | True | True |
| 8000000-numeric | row | 1 | 0.466589 | 0.523438 | inconclusive | 124.922 | False | False |
| 8000000-numeric | row | 5 | 2.21264 | 2.363592 | True | n/a | batch excluded | batch excluded |
| 8000000-numeric | column | 1 | 0.557658 | 0.623291 | inconclusive | 111.344 | False | False |
| 8000000-numeric | column | 5 | 2.675241 | 2.829042 | True | n/a | batch excluded | batch excluded |
| 8000000-logical | row | 1 | 0.535853 | 0.597087 | inconclusive | 124.391 | False | False |
| 8000000-logical | row | 5 | 2.32236 | 2.538987 | True | n/a | batch excluded | batch excluded |
| 8000000-logical | column | 1 | 0.621922 | 0.690783 | inconclusive | 123.641 | False | False |
| 8000000-logical | column | 5 | 2.793581 | 2.997908 | True | n/a | batch excluded | batch excluded |

Gate crossings: 4. No correctness-cost exception is automatically granted.

RSS is in bytes in JSON; sampled phase peaks can miss short peaks. Rprofmem records allocation traffic, not isolated live sorting workspace.
The full JSON retains every warmup, trial, profile, input mutation/output equality check, build/archive hash and threshold comparison.
Integrated workflow, backed, worker-memory/cache and consumer gates remain separate.
