# Permanent test inputs

These files are inputs to tests of current SCEPTRE behavior. They contain no
historical p-value, effect-estimate or random-stream output snapshots.
`provenance.json` records their SHA256 values; update it deliberately when an
input changes, alongside the reason and affected tests.

- `refactor/inputs.rds`, `cases.rds` and `backed-inputs.rds` describe synthetic
  counts, covariates and seven workflow configurations. The counts include a
  strong positive effect, a zero-count gene, failed-QC pairs, nonconsecutive
  retained cells and overlapping high-MOI guide assignments.
- `refactor/native-inputs.rds` contains explicitly constructed score arrays and
  draw lists, plus sampler settings. Stage labels describe the intended examples;
  tests also check independent effect/SE formulas and selected draw segments.
  Historical function-name strings in sampler metadata are unused: the current
  helper explicitly selects the production entry points.
- `aggregation-inputs-6a43fa5.rds` is a frozen set of 28 count-matrix/metadata
  examples introduced at 6a43fa5: duplicate triplets and metadata, stored zeros,
  logical/pattern matrices and rounding-sensitive counts. Keep these inputs
  frozen; add new examples in the test when extending coverage. Expectations
  are calculated from dense algebra during the current test, not saved outputs.
- `saved-objects/sceptre-0.99.0-21f9ea0.rds` is an in-memory object written by
  source-verified SCEPTRE 0.99.0 at 21f9ea098b69c444a884a49bc254de711968e3b5.
  It is the prepared low-MOI/complement/permutation toy analysis, before discovery.
  Its original bytes were preserved from the refactor characterization fixture.
  The current package must read it and resume equivalently to a newly prepared
  object. Do not regenerate it with current SCEPTRE: it is a versioned file-format
  example, much like a test input written by an older file-format implementation.
  No old executable or old expected result is required by this test.

To recreate the four mathematical/workflow input files in a fresh directory:

```sh
Rscript --vanilla tools/fixtures/generate-inputs.R /absolute/fresh/input-directory
```

This uses current test-input recipes, R, Matrix and withr, without loading SCEPTRE
or its original implementation. Compare the R objects with `identical(readRDS(a),
readRDS(b))`; compressed bytes can differ across R versions. Review input changes
rather than automatically replacing fixtures during CI. The bridge C++ file is
small test-only code for constructing and inspecting draw banks.
