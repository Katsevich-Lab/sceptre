# SingleCellExperiment interoperability

The public entry points are `import_from_sce()`, `set_grna_assignments()` and
`export_results_to_sce()`. Read `vignettes/single-cell-interoperability.Rmd`
for the raw-count, full-library, alignment and calling contracts. The optional
Bioconductor dependencies belong in an isolated library when qualifying this
candidate. This directory is source tooling, excluded from the package archive.

Run against an installed candidate package and the schema-3 pipeline companion:

```sh
R_LIBS=/absolute/isolated/library Rscript \
  tools/interoperability/sce-roundtrip.R \
  /absolute/pipeline-checkout /tmp/new-sce-roundtrip
```

The script uses all 100 responses and 1,000 cells in the bundled example. It
creates an SCE, imports raw counts, supplies explicitly labelled demonstration
guide calls, serializes a compatible pipeline bundle, and runs the pipeline
`stage.R` CLI in separate R processes for parameters, assignment, QC and discovery.
It checks unchanged full-library covariates and calls, exact direct/pipeline
result equality, original-SCE ownership, and saves returned SCE results plus
session/runtime identities. It does not run Nextflow, a container, a cluster,
or independently validate the scientific quality of the example caller.

The companion pipeline test `tests/external-assignments.R` checks preservation,
explicit reassignment, stale memberships and design invalidation. Both repositories
must be integrated together; a package-only install supports direct SCE work,
but an older pipeline recomputes assignments by default.

The routine regression file `tests/testthat/test-single-cell-interoperability.R`
checks raw-count ownership, library covariates after pair selection/QC, strict
metadata/cell/guide identities, invalid counts and calls, sparse calls,
vector-level assignments, high-MOI overlap, low-MOI exclusion, stale identities,
and owned results after serialization. It skips SCE cases when the optional
dependency is absent. Ordinary matrix import and external assignments retain
their existing required dependencies.

Implementation support is limited to in-memory raw-count assays. AnnData,
delayed/backed SCE import, upstream-normalization detection, and recovery
of upstream-filtered feature libraries are outside this slice.

## W05 terminal evidence

`validation.json` records the exact measured source/runtime and dependency
versions. The isolated install passed; all 56 targeted assertions passed without
failures or warnings. Pipeline preservation, explicit override, stale-call and
design-reset rejection passed with SCE absent. The fresh-process example retained
all 100 input response rows and 50 guide rows across 1,000 source cells, retained
986 cells after QC, returned two requested pair rows, and matched direct results
at tolerance zero. Optional-SCE absence checks passed.

The installed API source was `22087ac`; the subsequently reviewed source change
`a4431c8` spells one loop as `seq_len(2L)` to avoid a BiocCheck coding-practice
instance. The separately committed driver fix `9a01650` corrects subprocess
stdout/stderr logging; it does not change installed package code. An initial
driver run incorrectly interpreted captured stdout as a process status; its
replacement run completed with exit zero. All generated RDS/log artifacts remain
local under `/private/tmp/sceptre-w05-*`.

The final integrated archive/vignette build, matched BiocCheck category/instance
gate, scheduler/container qualification and independent scientific validation
remain program gates. These local API/serialization receipts do not supersede
them or qualify the foundation positivity issue reported by W02.

## Seurat companion

The separate `seurat-interoperability` vignette and `seurat-roundtrip.R` exercise
the optional SeuratObject 5 path. A single complete `counts` layer is required
per assay; split counts, partial layers and differing cell universes reject.
The adapter uses the existing external-call and pipeline contracts. Full Seurat
and SingleCellExperiment are unnecessary. Returned pair tables and provenance
live in `Misc(returned, slot = "sceptre")`. Initial SCE receipts above do not
qualify this newer companion; its versioned evidence is separate.

`seurat-validation.json` records the terminal new checks: 32 Seurat assertions
and 56 SCE regressions passed with no test failures or warnings. The actual
100-response/50-guide/1,000-cell example retained 986 cells and returned two
finite pair results through fresh pipeline CLI processes, matching direct
results at tolerance zero. A controlled fresh process with SeuratObject absent
passed ordinary import/calls and both optional-adapter guards.

The runtime was SeuratObject 5.4.0 on R 4.6.0/macOS ARM64, using rebuilt native
source with the accepted empirical cap and observed-score guard. The ambient
`sp` package emitted a built-under-R-4.6.1 warning; it is retained separately
from test warnings. An initial test-only data-free assay lacked its key and
was rejected by Seurat validity before the intended boundary. The corrected
fixture asserts validity first; the complete replacement run exited zero.
Original failed logs remain local. Final integrated archive/vignette, matched
BiocCheck, platform and scientific gates retain their separate authority.
