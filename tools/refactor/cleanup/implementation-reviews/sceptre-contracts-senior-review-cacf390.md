# Independent source-contract review

Reviewer: rf0_senior_review. Date: 2026-10-16.

Reviewed immutable commit `cacf3908d48e0b9352aa22813dfcee2f71f838ac` versus its parent. No production edits, build/test reruns or benchmarks. A lightweight independent parse confirms identical ASTs for all eight changed R files; the native numerical diff is comments only, and public plot Rd/source usages and aliases remain unchanged.

**No P0/P1/material P2 found.** Contracts were checked against actual helpers/callers for import/ID order and vector aggregation, backed parser delegation, raw/source/post-QC/NT-local coordinates, saved wrappers/56 slots, stored-position versus loaded-value counts, calibration default/selection order, sampler scope/RNG, directory/file effects, formula selection, numerical failure state and plot return direction. The comments document retained behavior rather than add scientific fallbacks or new backed capabilities. Existing full-input dense vector aggregation is accurately declared as an open gated implementation, not falsely described as already sparse.

One small wording correction is recommended at `R/aux_functions.R:5`: replace "at least 15 observed values" with "at least 15 distinct values". The code tests length(unique(column)), not number of rows/nonmissing observations. This is not a runtime blocker but improves the new source contract's precision.

AST parity establishes no runtime expression change, not that all prose is automatically correct or that source is pristine. This review does not replace final immutable installed/consumer/resource qualification for the combined branch.

## Commit-file SHA256

- `R/aux_functions.R`: `5079a9f9196d9050b438dcce3f5a48b4fac4e2e9bb3b90d6341d9cbefa35ba65`
- `R/getters_and_setters.R`: `32d34657856301e71212ef9211b96f54361a29f62f40b22a408dc93d479fea2e`
- `R/import_functs.R`: `bf9e43c4eb3a0cff6eda120101d071298554c5ba29e359121511890307c2e890`
- `R/neg_control_functions.R`: `f4fe402cc769cb8d823cf31f18a1b86c2c2a57f664a78aa20c19fce8cb6717dd`
- `R/pairwise_qc_functs.R`: `27f9ecb236ef0b56dbd2ebc880bc847201bd31b3d15fa1d91ad6eb516cc4a5e8`
- `R/preprocessing_functions.R`: `3a225440729bfe1fbee9ece14e88bd7a8c4fe6e38de873a786aae369da425792`
- `R/s4_class.R`: `3f66cd9c563e6e51ef1b0e838831b404bc52f08545494f48b11ea026f6539998`
- `R/s4_helpers.R`: `29a31262d8cae42b7dc48fd29d42897c58ba23a9a098be1c8db40a07a5d5d770`
- `man/plot-sceptre_object-method.Rd`: `a5e063660193ac29f55c15dabb9f0e894415c47607365bbae1032869e5faee86`
- `src/compute_nb_size_param_functions.cpp`: `c3cbea31f5ae745e30deda5e88cf4f5bb8d0a6c4bda890d28cb5b18b19d04146`
- `tools/refactor/STATUS.md`: `a170e2ba0518521fcf86c89c318fe00a4b564bb8cb53767bb03bf8effa5709ac`
