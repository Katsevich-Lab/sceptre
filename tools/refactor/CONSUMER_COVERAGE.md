# Consumer characterization

`pipeline-smoke.R` runs the pinned pipeline's original R scripts in fresh
processes, with an asserted installed sceptre path, real reattached `.odm`
files, isolated working directories, captured warnings and persisted observations.
The pipeline revision and all 11 script hashes must match `consumer-inventory.json`.
The harness never regenerates frozen package fixtures.

The expanded scope is 30 stages:

- High-MOI parameter defaults, guide pod planning, four threshold-assignment
  workers, assignment gather, full QC and trans-sentinel parameter handling.
- Low-MOI parameter defaults, guide planning, maximum-assignment gather, full QC
  and trans-sentinel parameter handling.
- Association preparation, two calibration workers and gather, one power worker
  and gather, and two discovery workers and gather.
- Nuclear cellwise QC, actual response-pod planning, two trial workers, two
  nontrial workers, and a nontrial shard containing only zero-passing pairs.

The nontrial zero-passing branch writes logical `NA` columns for `p_value` and
`log_2_fold_change`; the ordinary analyzed path writes numeric columns. This is
an explicit baseline schema distinction, preserved in comparison payloads.
Nontrial pair construction also retains the pinned script's guide-table-based
ordering and duplicate target groups. These are characterizations, not fixes.

```sh
R_LIBS=/path/to/selected-ondisc-library \
Rscript --vanilla tools/refactor/pipeline-smoke.R \
  /path/to/sceptre-library /path/to/pinned/sceptre-pipeline \
  tests/testthat/fixtures/refactor /new/output-directory
```

Run reference and candidate against the same dependency libraries and fixture
files. Check each report's statuses and warnings before comparing
`comparison.rds` with `all.equal(..., tolerance = 0)`. The payload includes
ordered result values/classes, facade flags, assignment memberships, cell views,
pair tables and cache contents. A report containing failed or blocked stages is
incomplete even if two incomplete payloads happen to compare equal.

The untouched sceptre source `21f9ea098b69c444a884a49bc254de711968e3b5`, installed
with source ondisc 1.3.5, passed all 30 stages on October 16, 2026. Its report is
`/private/tmp/sceptre-rf8-expanded-reference-final/report.json`. The installed
refactor candidate built from `f4aaab4` also passed all 30 stages, with an exactly
equal complete comparison payload, preserving decisions, schema and ordering.
Its report is `/private/tmp/sceptre-rf8-expanded-candidate-final/report.json`.

## Actual Nextflow execution

`pipeline-nextflow.R` archives the pinned pipeline revision without changing its
sources, prepares real backed inputs, and executes the complete standard and
nuclear workflows using the local executor. It retains archive/source hashes,
the complete runtime configuration, task logs and R-process package provenance.
It checks standard gathered results/global correction and nuclear raw shards,
then opens and collects the actual Arrow dataset across all shards.

```sh
R_LIBS=/path/to/selected-ondisc-library \
Rscript --vanilla tools/refactor/pipeline-nextflow.R \
  /path/to/sceptre-library /path/to/pinned/sceptre-pipeline \
  tests/testthat/fixtures/refactor /new/output-directory
```

This is an explicit local qualification environment: Nextflow 25.10.4, Java 17,
source ondisc 1.3.5, disabled containers/conda and job arrays, and two concurrent
tasks. The archived default `array = { params.array_size }` is rejected by this
Nextflow version; static arrays are unsupported by the local executor. The
external qualification config therefore overrides that exact selector with
`array=null`. This does not qualify the unmodified default array configuration,
historical containers, nf-test or hosted scheduler profiles.

The original nuclear pod planner calls `sample()` without a seed. Independent
unseeded workflow runs correctly produced different pod assignments and could
not be compared by shard filename. The qualification startup profile now fixes
seed 271828 in every R task, preserving a shared work plan before comparing
reference and candidate. No output sorting is used to conceal a changed plan.

Reference and candidate each execute all standard and nuclear tasks; paired
reports and payloads are under
`/private/tmp/sceptre-rf8-nextflow-{reference,candidate}-seeded/`. The response-to-pod
maps are identical and the complete ordered comparison payloads are exactly
equal (also passing the predeclared `1e-8` floating-point tolerance). Each lane
records 34 R-process provenance files for each of its two successful workflows;
the nuclear dataset readback contains all 32 emitted rows.

## Uncached book cells

`book-smoke.R LIB BOOK_CHECKOUT OUT` executes 15 unchanged R cells from the pinned
`import-data.qmd` and `ondisc.qmd` chapters in a fresh process with no knitr cache.
Constructor output files stay inside an isolated temporary directory. The
report retains exact code and line ranges, book/data hashes, warnings, package
paths and completed observations. It does not execute installation cells,
machine-specific cleanup commands or complete analysis chapters.

Reference and candidate each pass 14 cells with exactly equal completed
observations, including the Parse import and real ondisc construction/row
access. The Cell Ranger sparse-import cell fails identically in both: the
installed `sceptredata` example contains two `barcodes.tsv.gz` files, each only
20 bytes and zero decompressed lines. The existing importer finds these files
and `fread()` rejects their empty contents. Reports retain that failure and
return nonzero; these are not described as successful whole-book runs.

Paired reports are `/private/tmp/sceptre-rf8-book-{reference,candidate}-final/`.
Valid example barcodes, complete cis/trans/singleton cells, whole-book rendering,
external rd7 data, historical containers and hosted executor profiles remain
separate gates.
