# Screen and Analysis objects: first draft

A working first draft of slide 5 of the development plan, on `dev`. It is
meant to be reacted to: every choice below can change.

## What exists

| Part | In memory | On disk |
|---|---|---|
| Screen | `SingleCellExperiment`: response counts (`counts`) and cell covariates (`colData`); gRNA counts and each cell's top gRNA in the `gRNA` altExp; MOI in `metadata` | `screen.h5mu` (MuData), modalities `response` and `grna`; the gRNA table in `grna_targets.parquet` |
| AnalysisSpec | `settings`: pairs, formula, side, resampling, control group, thresholds, B1–B3, gRNA strategy and assignment settings, multiple testing, QC thresholds | tables in `tables/`, the rest in `state.h5ad` |
| StageRecords | one `StageRecord` per stage: `assignment`, `qc`, `design`, `pairs`, `calibration`; each with its `contents` and a hash of them | `state.h5ad`: per-cell values in `obs`, the design matrix in `obsm`, tables in `tables/`, everything else in `uns$state` |
| ResultTables | `calibration`, `power`, `discovery` | `results/NAME.parquet` |
| FitCache | per response: fitted coefficients and dispersion | `fits.parquet`, one row per response and term |

`R/analysis_objects.R` defines the classes and the lossless conversion to and
from today's `sceptre_object` (`as_analysis()`, `as_sceptre_object()`);
`analysis_layout` says where each slot lives. `R/analysis_files.R` writes and
reads the files (`write_analysis()`, `read_analysis()`). Nothing is exported
yet, and every analysis function still works on `sceptre_object`.

## What was checked

- The conversion returns an identical `sceptre_object` at every stage, from
  import to discovery.
- Writing and reading back returns an identical object at every stage of a
  low-MOI workflow (thresholding) and a high-MOI one (mixture assignment,
  calibration with `output_amount = 3`, power check), and resuming from any
  stage gives identical discovery results.
- Files cannot hold a formula's environment, data.table's internal pointer, or
  whether a table was a tibble (it comes back as a data.frame). None affects
  results.
- Python reads every file: the Screen with `mudata`, the state with `anndata`,
  the tables with `pandas`.

## Decisions to review

1. **Conversion, not replacement.** The analysis functions still run on
   `sceptre_object`; the Analysis wraps it. Moving each stage onto the Analysis
   is step 3 (modules).
2. **Stage boundaries.** Five stages, by which function sets the slots:
   assignment, QC, design (`set_analysis_parameters()`), pairs (pairwise QC)
   and calibration. Power and discovery produce only results.
3. **Hashes are recomputed, not stored.** Each stage's hash is of its own
   contents (`rlang::hash`), made when the Analysis is built, including after
   reading. Recording which inputs a stage used, so that a stale stage can be
   detected, is not done yet. `rlang::hash` is not promised to be stable across
   versions, so it should not be written to files as is.
4. **Where values go.** Per-cell values are `obs` columns (gRNAs per cell; the
   cells in use and with zero or two-plus gRNAs as logical columns); the design
   matrix is in `obsm`; tables are Parquet, which keeps their columns, types,
   missing values and row names; every other value is one JSON text in
   `uns$state`, written by `jsonlite::serializeJSON()` (17 significant digits),
   which keeps R types, names, order and missing values exactly. Formulas are
   stored as text.
5. **The settings JSON is R-shaped.** Python can read it, but each value carries
   its R type and attributes. If pysceptre reads the state, plain JSON would be
   friendlier.
6. **Scale.** The cell-index lists of gRNA assignments are in the JSON; at
   atlas scale they should move to Parquet or `obs`.
7. **Old pipeline fields kept.** `nuclear`, `nf_pipeline`, `integer_id`,
   `elements_to_analyze` and the never-set `mitochondrial_gene` are carried in
   `legacy`, to be dropped with ondisc.
8. **Dependencies.** `SingleCellExperiment`, `SummarizedExperiment`,
   `S4Vectors` and `rlang` in Imports; `MuData`, `MultiAssayExperiment`,
   `anndataR`, `arrow` and `jsonlite` in Suggests, needed only for files.

## Gaps found upstream

- **MuData 1.16** reads a missing value in a text column back as `"NA"`, and
  in a numeric column as `NaN`. `read_screen()` works around both for cell
  covariates; the fix belongs in MuData. A draft issue is in the overnight
  summary.
- **anndataR 1.2** returns vectors as one-dimensional arrays, sorts list
  entries by name, and drops values it cannot write with only a warning.
  Storing settings as one JSON text avoids the first two; warnings are turned
  into errors.
