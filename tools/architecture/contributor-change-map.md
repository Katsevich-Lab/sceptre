# Find the right place for a change

Begin with a concrete behavior, identify the owner and function boundary, and include the tests that establish the intended contract. An extension point is a place to start reviewing a change, not a promise that a new method needs no coordinated edits.

| If you are changing | Start here | Contracts to check |
|---|---|---|
| Imported counts or metadata | [Import](../../R/import_functs.R), [storage](../../R/storage.R) | Aligned cells, canonical IDs, count semantics, caller ownership, content identity and invalidation. |
| Scientific choices or requested hypotheses | [Architecture](../../ARCHITECTURE.md#boundary-contracts), [AnalysisSpec fields](../components/schema3/OBJECT_FIELDS.md#analysisspec) | Canonical requests, validated choices, preparation dependencies and result provenance. |
| Assignment or QC | [Preparation lifecycle](../components/schema3/ARCHITECTURE.md#preparation-lifecycle-and-fitting) | Source versus analysis positions, stage identity, retained fitting cells and independent cache invalidation. |
| Response fitting | [Fitting](../../R/fitting.R), [score preparation](../../R/score_preparation.R) | Compact coefficients and dispersion, actual fitting scope, method identity, cold/warm behavior and scientific validation. |
| CRT perturbation fits | [Perturbation fits](../../R/perturbation_fits.R) | Exact probability reconstruction, actual treatment/design/cell identity, fresh draws and calibration bypass. |
| A statistic or tail calculation | [Workers](../../R/association_workers.R), [native entry](../../src/low_level_full_test.cpp), [tails](../../src/inference_tail.cpp) | Effect/statistic/output shape, numerical boundary behavior, RNG and draw ordering, complete-family adjustment. |
| Result handling | [Result processing](../../R/result_processing.R), [result records](../../R/result_records.R) | Row identity, failed-QC rows, global adjustment, owned snapshots, provenance and currentness. |
| Worker execution | [Execution boundary](../../ARCHITECTURE.md#concrete-extension-seams), [distributed value protocol](../../R/distributed_execution.R) | Complete-family scope, errors versus NULL replies, process-local resources and atomic collection. |
| Saved-object compatibility | [Legacy field decisions](../components/schema3/FIELD_DECISIONS.md), [migration](../../R/schema_migration.R) | Physical saved layouts, historical results, unknown provenance and caller preservation. |

## A bounded contribution sequence

1. State the expected public behavior and supported scope.
2. Find the authoritative owner and existing function boundary.
3. Read the associated contract and regression files in [the architecture table](../../ARCHITECTURE.md#boundary-contracts).
4. Make the smallest coordinated implementation, tests, and documentation change.
5. Run the relevant checks using [the contribution instructions](../../CONTRIBUTING.md) and [maintainer check guide](../checks/README.md).
6. Describe behavior, validation, and remaining limitations for review.

Changes to model fitting, calibration, power, or numerical results need scientific validation appropriate to that change. Passing unit tests alone does not establish a new method's statistical validity.

## Setup and checks

The [package README](../../README.md) gives installation and optional-backed-dependency instructions. [CONTRIBUTING](../../CONTRIBUTING.md) gives the formatter and package-check commands. Optional ondisc paths have an explicit [dependency contract](../hardening/ONDISC_CONTRACT.md).

The existing setup instructions are rendered here for review. A fresh-checkout contributor exercise and an executed worked extension remain to be verified before contributor readiness is signed off.
