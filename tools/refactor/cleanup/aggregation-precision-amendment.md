# B6 platform-precision amendment

Hosted Linux CI at 0bb30e8 exposed a real replacement defect. The original
21f9ea0 function first calls `as.matrix()` and then `Matrix::colSums()` on dense
selected rows. On Linux this uses extended accumulation; the sparse replacement
uses doubles. `(2^53,1,1)` and `(.1,.2,.3)` distinguish the paths. The former test
also incorrectly froze the Mac result as universal. Neither a tolerance nor a
platform skip resolves the production difference.

The initial 88ed37c source review and complete Mac exact/resource evidence are
preserved unchanged, but their cross-platform acceptance is superseded. The new
private helper uses the same dense Matrix sum on complete selected guide rows,
in column blocks targeting 2^20 elements (8 MiB for doubles). One entire column
is the minimum: more than 2^20 selected rows can exceed that workspace target.
Dense vector-by-cell output/transpose and sparse row-selection workspace are
additional. Splitting rows and adding partial sums is prohibited because it
changes accumulation order/precision. No custom compensated/native sum is added.

The dplyr selector, row repetition/order, triplet normalization, logical FALSE
values, first-seen vectors, cell names and t(vapply()) assembly remain. Tests
execute an unchanged test-only copy of the original function on this R/Matrix
platform. Its bytes and the 28 previously captured inputs are hash-pinned in
separate fixtures, outside the immutable RF0 goldens. Another 200 seeded triplet
cases and a 1025-by-1025 block-boundary case compare exactly; zero/one columns
are explicit. No RF0 fixture or tolerance is changed.

Before acceptance: independently review the implementation and reference
provenance, run these tests in fresh installed candidate namespaces including
hosted Linux, and repeat the previously declared bounded B6 resource protocol
against matched fresh original/candidate release builds. Preserve every earlier
trial and any new failure. The gates/inputs/repetitions/safety limits remain
those in vector-aggregation-resource-plan.md; the changed summation path has no
resource exception. The separate accepted numeric-triplet ownership cost does
not apply here.

## First repaired-source resource attempt

At source 9631ef2 the fresh installed original/candidate replay passed all 228
outputs, and focused installed tests passed 595 checks in 17 blocks with no
failure/error/warning/skip. Independent source review closed. The first renewed
resource attempt is preserved at
`/private/tmp/sceptre-b6-precision-9631ef2-resource-run1` with its complete log.
The large numeric setting passed timing/RSS (ratio 0.913; original/candidate OS
peak 1426.8125/664.765625MiB). The medium numeric triplet setting passed RSS but
crossed the 10% timing gate: 20-call medians 2.086/2.399s, ratio 1.150. The harness
stopped as declared; logical/actual-data settings were not reached. No exception
was accepted, and none of these trials are discarded or relabeled.

Source inspection identified avoidable overhead: even a whole selected guide
group already within the 2^20-element target allocated another output vector,
column sequence and sparse column subset. Such groups now directly execute the
same dense `Matrix::colSums(as.matrix(counts))`; only larger groups use column
blocks. Dimensions multiply in double. Complete row order/precision and the
workspace target are unchanged. This changed source requires fresh installed
and resource qualification; a repeated run is justified by the implementation
change, not by selecting favorable timings.
