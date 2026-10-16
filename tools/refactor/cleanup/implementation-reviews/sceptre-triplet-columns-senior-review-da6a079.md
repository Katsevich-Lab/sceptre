# Static triplet-column ownership review

Reviewer: rf0_senior_review. Date: 2026-10-16.

**No introduced P0/P1/material P2 in production copy reduction `da6a079dc86bc18d7dd54e4ba590820eacd921ca`.** Compared its four-file unit with its parent; no native/binding changes. This closes static column-ownership review, not resources or all retained dense-conversion correctness. The separately reproduced baseline symmetric-dense defect is recorded in `/private/tmp/sceptre-symmetric-dense-disposition.md` and does not arise from these ownership flags.

The three flags follow constructors rather than NAMED/refcount heuristics: supplied sparse slots remain borrowed and get data.table::copy; generated pointer coordinates, replacement logical-one values and freshly coerced dense sparse columns may be sorted directly. The local list/setDT owns all sorted columns. No numeric arithmetic is used to copy, so negative-zero/NaN payload bits are preserved. No scientific values, stored logical-FALSE presence rule, duplicate ordering, dimensions/dimnames, sparse output orientation or pointer algorithm changes. No speculative fallback or version branch is introduced.

The installed Matrix1.7.5 public matrix→Tsparse method and .m2sparse.checking/.m2sparse route were inspected. The retained versioned C source hash67805fdd88ccdf06a4c1d888d55df3dec92031d27b20d5410a23cec39ab657b3 allocates fresh i/j/x for triplet storage (1969/1975/2019), even when the constructor chooses a structured shape. That freshness proves ownership only; structured symmetry semantics require the separate correction. Stale rationale line links were corrected in docs-only f3d741a without changing production da6a079.

New tests cover both output orientations, compressed/dense/logical input and caller aliases, explicit FALSE semantics, and bytewise -0/nondefault-NaN preservation. The worker reports1095 focused installed expectations on a fresh immutable release build; no tests/builds were rerun here. The profiler now traces explicit copy as well as table construction/sorting, without changing gates or primary timings.

Initial8841efa resource results remain immutable, include memory crossings and are not a pass. New source must receive its own matched release resource run after the separate dense correction's scope is settled. This review grants no performance gain, threshold waiver or ownership-cost exception. Sparse/backed final resources and general combined-source qualification remain separate.

## File SHA256

- `R/storage.R` at da6a079: `3755dea9c81629b0140b03230a09d6359e57169ea02580ce4f19c1da9b4df06e`
- `tests/testthat/test-refactor-storage-ownership.R` at da6a079: `59ddd65bb8735bc3a2d9c3d99f300f1b5b35aba67429642642808b859d7041e9`
- `tools/refactor/TRIPLET_RESOURCES.md` at da6a079: `9e4947ce4a53770ba3b2db379c30566fd548ede6dd7622529ef1f70b031dae76`
- `tools/refactor/triplet-resource-worker.R` at da6a079: `95bed39fc9404fb163590668f9d460fcdf018209e7c36b36658c4f50d355719e`
