# Dense general-expansion source review

Reviewed exact commit `bc53bc1b09abff72e5d4745c46d173bc1f77ec17` against its parent, read-only. No R process, build, test or benchmark was run for this review.

**Outcome: no P0/P1/material P2 source finding in this bounded repair.** It fixes the previously reproduced accepted-input loss by forcing general structure before sparse coercion, rather than expanding a triangle after auto-symmetry has already discarded a distinct lower value. It uses the supported public Matrix coercion route; no bespoke slot/validity bypass or private Matrix dependency is introduced.

The coercion is confined to the base-matrix branch. Sparse and ondisc branches and their static ownership decisions are unchanged. The resulting sparse coordinates/values are fresh under the already inspected Matrix allocation implementation, and the existing owned sort remains appropriate. Documentation candidly identifies an additional dense-value copy and makes no resource-pass claim.

The tests use independent 2x2 symmetric and slightly asymmetric full matrices with matching dimension labels, the condition that exposed auto-symmetry. They verify both orientations, exact tiny off-diagonal values, slot order/classes and serialized caller/alias nonmutation. The public import checks exercise getters, cell totals, nonzero counts and mitochondrial fractions against full independent matrix calculations; they do not turn the buggy baseline result into an expectation.

Runtime installed-source qualification and dense-copy cost/resource measurements remain open. This is source closure of the narrow repair, not acceptance of the complete storage resource unit, arbitrary Matrix versions, or the final integrated package.

## Source SHA256

- `R/storage.R`: `f5d33db927292a1f562f2128c5fd37e699b6be93fc4c43a952acefcafa7f1398`
- `tests/testthat/test-refactor-storage-dense.R`: `f6f46a0da0c1c53fa8ac618d77a4a417b77fd5faba54a76e8245c2443be22031`
