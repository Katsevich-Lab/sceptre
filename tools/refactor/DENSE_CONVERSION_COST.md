# Dense general-structure correction: bounded cost characterization

Dense symmetric and tolerance-near-symmetric inputs must retain every actual
entry before conversion into a general sparse matrix. The supported public
Matrix route is dense matrix to `generalMatrix`, then `TsparseMatrix`. It
requires an additional dense payload copy; sparse/backed inputs do not use this
route and are never densified. A hand-built wrapper is deliberately avoided.

Measure this distinct correctness cost against original
`21f9ea098b69c444a884a49bc254de711968e3b5` and the accepted repaired source pin,
with fresh matched release `-O2` libraries, source/archive/package/DLL hashes,
source ondisc 1.3.5 and all four thread controls plus data.table threads at one.
The small symmetric/near-symmetric tests separately demonstrate the corrected
behavior. Performance cases are rectangular, so their sparse outputs must be
exactly identical in reference/candidate; do not compare defective triangular
baseline results as if they were an equivalent mathematical operation.

Use two immutable numeric dense fixtures, 1000 by 1500 and 2000 by 3000, with
10% nonzero entries at R linear positions 1, 11, 21, ... . Nonzero values are
`((seq_along(positions) %% 97) + 1) / 7`; row/cell names are deterministic.
Payloads are 12,000,000 and 48,000,000 bytes. No random generation or scientific
golden changes occur. Both row and column outputs are measured separately.

Use the same conversion clock, fresh independent input loading, pre-call GC,
input/output live markers, 60 ms untimed RSS padding, exact sparse slot/dimname
comparison, OS whole-child peak and timestamped 20 ms phase sampler as the
triplet protocol. Run eight untimed correctness conversions/four comparisons,
eight warmups, four measured process pairs per setting (32 children) with
original-first/candidate-first/candidate-first/original-first ordering, and eight
additional allocation-profile children. Retain all trials and Rprofmem/heap
records. Report measured dense-copy allocation and RSS deltas per input element
as well as per stored sparse entry; allocation traffic is not peak live memory.

Use the same 8 GiB headroom check, 1536M R vector heap and conservative 2 GiB
RSS stop. One child runs at a time in a confirmed quiet window. A safety stop or
missing case remains incomplete. No timing-only batches are planned for this
bounded characterization: medians below one second are explicitly inconclusive
as ratios; retain absolute timings. Report the existing 10%/64 MiB RSS and 10%
meaningful-time comparisons transparently, with no automatic cost waiver or
claim that this replaces the independent full triplet/workflow/backed gates.

Run the existing `triplet-resource.py` phases with `--recipe dense` and the
same required library/build/archive arguments. Its recipe/manifest distinguishes
this bounded dense characterization from the approved triplet experiment.
