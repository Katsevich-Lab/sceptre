# SCEPTRE refactor handoff

Download [the readable tour ZIP](sceptre-refactor-tour.zip), extract it and open
`walkthrough.html`; its links open the bug appendix. No render tool is needed to
read it. The [tour source](walkthrough.qmd), [appendix source](bug-fixes.qmd) and
lasting [architecture guide](../../ARCHITECTURE.md) are also available.
The public API and
56-slot saved-object format are preserved. The refactor adds no new inference
method, storage backend or replacement QC policy. Normal numerical behavior is
checked against original commit `21f9ea098b69c444a884a49bc254de711968e3b5`.
Reviewed bug domains and improved worker-error reporting are explicit exceptions.
From a repository clone, `gh pr checkout 234` selects the review branch. The
[validation receipt](validation/README.md) gives the final measured outcomes and
source pins. The portable tour is a ready-to-read entry point for Tim.

## Material to keep

Permanent package tests check current mathematical behavior, workflow results,
sampling probabilities, warnings, indexing, cache reuse and supported execution
modes. The older saved object is a lasting file-format input, not an old-output
snapshot. Input provenance and maintenance instructions are in
[tests/testthat/fixtures](../../tests/testthat/fixtures/README.md).
Current-source native safety checks live in [tools/checks](../checks/README.md).
The maintained diagram lives in [tools/architecture](../architecture/README.md).

This directory alone is temporary. It retains the tour, exact manual comparison
suite and the explicit review commands below for maintainer review. Routine CI
never calls it. All developer directories are excluded from `R CMD build`.

## Full historical evidence

The complete pre-handoff tree, original specification, plans, frozen reports,
raw evidence and failed attempts are preserved at the published tag
[`refactor-evidence-2026-10-16-e1f8109`](https://github.com/Katsevich-Lab/sceptre/tree/refactor-evidence-2026-10-16-e1f8109/tools/refactor).
The [frozen acceptance report](https://github.com/Katsevich-Lab/sceptre/blob/refactor-evidence-2026-10-16-e1f8109/tools/refactor/FINAL_ACCEPTANCE.md)
records measured source pins, actual Nextflow/book scope and accepted resource
costs. The [specification](https://github.com/Katsevich-Lab/sceptre/blob/refactor-evidence-2026-10-16-e1f8109/tools/refactor/REFACTOR_SPEC.md)
remains available there. Current validation belongs to the handoff receipt;
historical numerical/resource measurements retain their original commits.

787 execution/planning/evidence files were removed from the active tree after
publishing that archive. This reduces the review diff; their Git history remains.
After maintainer acceptance, archive the final handoff tree before retiring this
directory. Keep all permanent tests, input recipes, native checks and architecture
material. Remove this directory's links from contributor guidance at retirement.

## Compare original and current SCEPTRE

Use a committed checkout, a C++ toolchain and the package dependencies. Install
`pkgbuild`, `pkgload`, `testthat`, `digest` and `jsonlite` in your developer library.
For backed qualification, follow the permanent README's exact guarded ondisc
installation; set `R_LIBS` to that library and require it explicitly:

```sh
R_LIBS=/absolute/dependency-library SCEPTRE_REQUIRE_ONDISC=true \
  Rscript --vanilla tools/refactor/compare.R /absolute/fresh/comparison-output
```

Without ondisc, in-memory comparisons run and backed tests skip; the report
states that scope. An installed unqualified ondisc is rejected before running its
kernels. An optional second argument selects a candidate Git ref; HEAD is the
default and uncommitted changes are excluded. Output must be a fresh directory
outside this checkout. The command archives both refs, freshly installs the
candidate and runs the verified original in a separate process. Read
`comparison.json`, `checks.csv` and the retained logs. Failures return nonzero and
record a failed status, stage and error. No candidate output regenerates an
expected result. This command does not qualify performance or actual Nextflow.

## Read or reproduce the illustrated tour

The handoff includes [a portable tour ZIP](sceptre-refactor-tour.zip): download,
extract and open `walkthrough.html`; keep `bug-fixes.html` beside it. Both pages
embed their images/scripts and link source excerpts to the recorded Git commit.
The ZIP's `READ-ME.txt` and `render.json` identify that source and render tool.
To reproduce it, install Quarto, R packages `knitr` and `jsonlite`, and the `zip`
command, then run:

```sh
Rscript --vanilla tools/refactor/render.R /absolute/fresh/tour-output
```

The optional second argument again selects a committed Git ref. Rendering reads
snippets from an immutable archive; it does not load an installed SCEPTRE or run
the comparison. Authoring directly with `quarto render walkthrough.qmd` remains
possible, but that renders the working checkout. Editable diagram sources are
Graphviz `.dot` files. Rendered working HTML is ignored; the compact handout ZIP
is retained for review and is excluded from the package.

## Adversarial test check

The explicit command below checks genuine current outputs first, then verifies
that nine deliberately wrong outputs are rejected by their intended assertions.
It exercises the same helpers used in permanent tests and does not alter the
package or any expected fixture:

```sh
Rscript --vanilla tools/refactor/check-test-strength.R /absolute/installed-library /absolute/fresh/check-output
```

## Follow-up ownership and acceptance

Eugene, as Lead Maintainer, owns triage of the bounded accepted costs in
[#235](https://github.com/Katsevich-Lab/sceptre/issues/235) and dependency-release
coordination in [#236](https://github.com/Katsevich-Lab/sceptre/issues/236).
Neither acceptance authorizes larger resource costs or enables backed CRT.
An inherited numerical edge case discovered during handoff is separately tracked
in [#237](https://github.com/Katsevich-Lab/sceptre/issues/237): fully tied two-sided
empirical tails can exceed one. Its original statistical formula is preserved in
this refactor; the new probability-bound tests concern the nondegenerate toy
workflow inputs. Eugene owns scientific triage of that separate correction.
Existing scientific and pipeline roadmap work stays separate: mixture replacement
[#220](https://github.com/Katsevich-Lab/sceptre/issues/220), pairwise QC
[#225](https://github.com/Katsevich-Lab/sceptre/issues/225), and pipeline plotting
[timothy-barry/sceptre-pipeline#6](https://github.com/timothy-barry/sceptre-pipeline/issues/6).
The final handoff also freshly executes all 30 original and candidate pipeline
R-script stages; their ordered payloads match exactly. This checks actual
standard/nuclear phase scripts, distinct from the earlier full Nextflow run.
Inherited nuclear-script options and the book example are owned for triage in
[#238](https://github.com/Katsevich-Lab/sceptre/issues/238).
The recorded book qualification is selected cells (14/15, with the inherited
empty-barcode failure); pipeline qualification uses matched local arrival/seed
settings. Full-book, container/array/hosted-scheduler validation belongs to their
maintainers' next qualification work, not a new feature promise from this PR.

Maintainer review of the final diff precedes merge/release. No one has been
messaged, assigned a review, or represented as having approved the refactor by
this handoff. This is the reproducible, reviewable delivery prepared for that step.
