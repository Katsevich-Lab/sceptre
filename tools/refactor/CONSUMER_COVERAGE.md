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

The latest completed paired qualification in this guide is original
`21f9ea098b69c444a884a49bc254de711968e3b5` versus candidate
`28df3907f43311efc70908d50b7636c7e6e89db8`, using verified fresh matched-O2
installs and source ondisc 1.3.5. Both pass all 30 stages with exactly equal
complete ordered comparison payloads. This is an earlier source checkpoint,
not qualification of later cleanup. [Complete evidence](results/consumer-qualified-28df390/README.md).

## Actual Nextflow execution

`pipeline-nextflow.R` archives the pinned pipeline revision without changing its
sources and runs complete standard/nuclear workflows using real backed inputs,
an isolated local executor and fresh task startup profiles. It records source,
configuration, logs, installed package paths and result payloads, then reads the
actual Arrow dataset across every emitted nuclear shard. The four arguments are
package library, pinned pipeline checkout, frozen fixture directory and new
output directory.

The paired 28df390 study retains three distinct runtime controls:

| Configuration | Execution | Ordered comparison |
| --- | --- | --- |
| Two concurrent tasks, fixed seed 271828 | Both workflows complete for both builds | Fails: baseline factor levels/tied rows depend on worker arrival; all five raw worker payloads are exact |
| queueSize=1/maxForks=1, fixed seed | Both workflows complete | Fails: serialization alone does not impose input emission order |
| Same serialized plan plus `fair=true` | Both workflows complete | Passes exactly, including identities/classes/order and global correction |

No published result is sorted or modified to conceal either failed control.
For exact final qualification use the unchanged arrival-controlled runner in
`results/consumer-qualified-28df390/runners/pipeline-nextflow-arrival-controlled.R`,
with its adjacent `consumer-inventory.json`, and fresh output directories. It
makes the complete source pin, package path, work plan and emitted data
reviewable. The generic two-task runner remains a separate control.

The qualified runtime is Nextflow 25.10.4/Java 17, source ondisc 1.3.5, one CPU/task,
containers/conda disabled, cached offline plugin and `array=null` for the exact
archived process selector. The default dynamic array closure is rejected by
this Nextflow version and local execution does not support job arrays. Each R
task starts with seed 271828 because the original response-pod planner samples
without a seed. Actual paired pod maps and arrival plans must compare equal;
independently unseeded runs are not comparable by shard name.

This qualifies the explicit local plan. It does not qualify uncontrolled gather
order, unmodified arrays, historical containers, nf-test or hosted scheduler
profiles. [All successful and failed controls](results/consumer-qualified-28df390/qualification-comparison.json)
remain tied to their original source. Final cleanup requires a fresh paired run.

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

The matched-O2 28df390 reports and exact payloads are retained in
`results/consumer-qualified-28df390/`; these supersede the earlier build claims
without claiming the known failed cell passed.
Valid example barcodes, complete cis/trans/singleton cells, whole-book rendering,
external rd7 data, historical containers and hosted executor profiles remain
separate gates.
