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
`/private/tmp/sceptre-rf8-expanded-reference-final/report.json`. Refactor-candidate
comparison remains a separate run, not implied by this reference result.

Full Nextflow execution, historical containers and uncached book workflows are
separate gates. These fresh-process script results alone do not close them.
