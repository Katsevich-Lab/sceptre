# Independent interpretation of final triplet diagnostic

**The diagnostic establishes extra preparation allocation traffic from protected triplet conversion. It does not establish the cause of the final actual Papalexi cold-serial +162.30 MiB RSS crossing.** No new source-correctness finding or architectural redesign is compelled by these results.

Reviewed final diagnostic `comparison.json` SHA256 `96e0e0e6c62d4565723df0091170567fb14fa048b4031584365c446fe2f092f3`, its reviewed driver, raw allocation profiles and actual targeted resource comparison. Candidate pin is `0a9f9d879db6a52e3daef237daa1711b9d542976`; actual reference is `21f9ea098b69c444a884a49bc254de711968e3b5`. I independently recounted profile categories, verified raw profile hashes, recomputed paired RSS deltas, and checked scientific/input/configuration anchors and before/after dependency snapshots. No R, profile, benchmark or source mutation was run by this reviewer.

## What the RSS contrast shows

One warm-up per role is excluded; three fresh measured pairs have 17–18 analysis RSS samples per process and identical scientific hashes/configuration.

| Analysis-only sampled main RSS | MiB |
| --- | ---: |
| Diagnostic production maximum | 570.52 |
| Unsafe borrowed-path maximum | 582.80 |
| Paired production-minus-borrowed differences | -15.30, +2.55, +6.14 |
| Median paired difference | +2.55 |
| Actual acceptance reference maximum | 504.97 |
| Actual acceptance candidate maximum | 667.27 |
| Actual candidate-minus-reference crossing | **+162.30** |

The diagnostic has no consistent large RSS reduction from borrowing. Its exact-production-source control itself peaks 96.75 MiB below the actual candidate experiment, so absolute acceptance peaks are not reproduced by these controls. Both arms deliberately source/byte-compile a replacement function, which changes their execution/heap context relative to the installed acceptance run. Three pairs, sampled lower-bound peaks and heap history do not identify the missing cause or prove zero copy contribution. The unsafe arm also changes sharing, caller mutation and sorting behavior: it is neither the original package nor an acceptable candidate remedy.

## What the allocation profiles show

The separate preparation profiles record 123,225,368 versus 73,534,952 allocated bytes: **49,690,416 bytes (47.39 MiB) extra traffic**, exactly 12 additional large vector allocations. Their sizes/counts match two sets of i/j/x copies for the response triplet (1,390,126 stored entries) and guide triplet (162,681), including allocator headers/alignment. Source tracing confirms two conversion contexts: column access for cell covariates, then row access for retained stores. The storage-stack delta is 49,706,792 bytes (47.40 MiB); overlapping stack categories and other small allocations account for the slight difference. This is allocation traffic, not retained/live bytes, process RSS or an unavoidable minimum.

The `copy_stack` lexical category is **misleading for copy attribution**: it reports the same 15 events/35,736 bytes in both arms, but raw profiles record the large changed allocations as `"<Anonymous>" "set_matrix_accessibility" ...`, without the literal `"copy"` or `"data.table::copy"`. The named-copy matches are small table operations in assignment/QC. Equality of this category does not mean the triplet copies allocate nothing; the raw large-vector contrast establishes the opposite.

All 56 per-slot `object.size` values are identical (summed 22,534,280 bytes in each arm; summation can count shared data repeatedly). At the preserved post-drop/pre-analysis GC point, live Vcells are approximately 45.3 MiB in both arms, differing by only 31 cells; Ncells differ by one. Vector heap trigger/high-water values differ (118.0/110.5 versus 94.4/78.2 MiB). This supports transient allocation/heap-history effects rather than a larger retained facade. These metadata come from separate profiled processes and do not explain or measure the actual acceptance RSS difference, native allocations or all sharing.

## Bounded disposition

Keep the nonmutation/ownership repair. Record its demonstrated +47.39 MiB preparation allocation traffic as a bounded correctness cost, and retain the **+162.30 MiB actual workflow RSS gate failure as a separate, causally unexplained measurement**. Do not relabel it as an established copy minimum, a safe borrowing opportunity or automatically covered by the earlier 909 numeric-triplet exception.

The remaining decision concerns acceptance of this measured workflow memory overrun under an explicit, scoped disposition versus requiring the original RSS gate. It is not a new numerical/API/ownership design decision. This diagnostic provides no reason to open an optimization or additional review program on its own, and grants no exception. Final acceptance and any exception remain with the task owner/user's existing authority.

Signed `/root/whole_package_review`, 2026-10-16. Companion `/private/tmp/sceptre-final-triplet-diagnostic-interpretation.json` retains exact hashes, independent arithmetic and anchor checks. Conclusions apply to this workload/platform and these measured settings only.
