# B6 bounded fast-path source review

Exact commit `783a2501295acc0b3fddab7bd2a615e545e5ae37` on top of the independently reviewed 9631ef2 precision repair. Reviewed the five production lines and qualification amendment, read-only; no R process/build/test/benchmark.

**No P0/P1/material P2 source finding.** The direct branch computes the extent in double and admits only selected groups already within the 2^20-element dense target. It calls precisely the same as.matrix then dense Matrix::colSums as each generic column block, preserving complete selected-row order, repeated metadata and platform accumulation. It avoids the additional sparse column subset, result allocation and assignment for the ordinary small-group case; no row-partial sums, alternate accumulator or statistical rule was added. Large-group block code is unchanged.

Zero/one row or column cases naturally take the direct branch with the original dense dispatch/shape handling. One minor nonblocking cleanup residue remains: after the direct branch, the old zero-column early return is unreachable for valid dimensions. This is not a correctness or resource finding and does not invalidate the running immutable build.

The completed 963 resource attempt retains its numeric-triplet timing crossing (2.399 vs 2.086 s, ratio 1.150) and successful primary/RSS observations. Those old observations do not qualify 783. Root's source-evaluated checks are not a fresh installed/hosted/resource acceptance claim. New immutable installed, same-environment exact, hosted Linux and renewed bounded resource gates remain required; no triplet-cost waiver applies to B6.

Source SHA256:

- `R/import_functs.R`: `aeeabf3386c4e65ea65d99cba319ab4db97394c688110e99ccc798534c568372`
- `tools/refactor/cleanup/aggregation-precision-amendment.md`: `5a1590c1440ef82c0a303b408b1debe0d3a53fb4105c77023bc33137fd71c3c0`
