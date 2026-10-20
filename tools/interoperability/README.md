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

Implementation support is limited to in-memory raw-count SCE assays. AnnData,
Seurat, delayed/backed SCE import, upstream-normalization detection, and recovery
of upstream-filtered feature libraries are outside this slice.
