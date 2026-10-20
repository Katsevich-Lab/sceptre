# Independent scientific validation

This harness generates **actual outcome counts and treatment assignments** and
calls the installed package's public import, parameter, assignment, QC and test
workflow. It does not synthesize Uniform p-values, compare against frozen p-values,
or treat direct/distributed agreement as statistical validation. The public method
contract is the [book's methods chapter](https://timothy-barry.github.io/sceptre-book/methods.html)
and [calibration chapter](https://timothy-barry.github.io/sceptre-book/run-calibration-check.html).
No new method or resampling budget is selected by this harness.

## Run

From the source checkout root, install the exact revision into a **new isolated
library**, with dependencies read from existing libraries. In the overnight
program acquire `work/overnight-2026-10-19/resources/local-heavy` via atomic mkdir
before compilation or simulation; use the absolute shared program path, write
owner/PID/command/start, release on exit, and cap BLAS/OpenMP/MAKE at two CPUs.
Do not install into ordinary user libraries or reuse another worker's mutable
installation. Example commands after the lease is acquired:

```sh
mkdir -p /private/tmp/w03-example-library
MAKEFLAGS=-j2 R CMD INSTALL --library=/private/tmp/w03-example-library .
Rscript tools/validation/check-harness.R
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
 R_LIBS=/private/tmp/w03-example-library Rscript tools/validation/run-scientific.R \
 /private/tmp/w03-pilot 40 600 all
```

Arguments are new output directory, positive replicate count, cell count (default
600), comma-separated case IDs (default all), optional candidate setup R file (or `-`),
and optional first replicate index (default1). Non-overlapping bounded batches
can use first indices1,201,401,601,801 with200 replicates each. Seeds and recorded
replicate IDs use those global indices; each batch remains independently terminal.
Never count the repeated first40 pilot seeds again when reporting the1000-dataset
combined estimate. `summarize-batches.R` rejects duplicate case/replicate/endpoint
identities and mismatched scientific source/binary/settings before aggregation.
Existing output directories are rejected. A setup file can set the option
`sceptre.validation.configure` to a function of `(object, case)` receiving the fully prepared object after QC and returning an
object with an explicitly selected experimental method. The script records the
setup hash; candidate APIs/unsupported cases must be documented with that setup.
For a separate experimental entry point, set
`sceptre.validation.run` to a function of `(object, case)` returning its actual
result table with `response_id`, `pass_qc`, `p_value`, and optional `stage` columns.
The harness controls the test seed. This hook invokes a real alternative test;
it must not synthesize or modify returned p-values or change the truth generator.
The setup is copied into the output directory for reproducibility.

`check-harness.R` verifies generator replay/change of seeds, low/high MOI structure,
truth values, and aggregation of errors, QC exclusions and out-of-range p-values.
The routine `test-scientific-sentinels.R` independently generates a strong positive
and negative count effect and checks public default inference for detection and
direction. It is a regression sentinel, **not** a type-I error certification.

## Truth and scope

Each case/replicate regenerates 600 independent cells (configurable). `z ~ N(0,1)`
and `batch ~ Bernoulli(.5)`. Treatment is Bernoulli with log odds -1.4, or
-1.4 + 1.3*z + .6*batch in confounded cases. Low MOI assigns exactly one guide:
treatment A, otherwise B with probability .3, otherwise one of four NT guides
uniformly. High MOI adds five independent Bernoulli(.3) guide indicators, so
multiple guides overlap. Present guide counts are 50; absent counts are zero;
thresholding at 5 recovers exact assignment. Assignment error/mixture fitting,
multiple targeting guides, correlated cells/donors and real screens are outside
this suite's validation scope.

Three independent response draws conditional on the same cells have NB size 3,
log mean `1 + .8*z + .4*batch + beta*X`, for beta 0, log(.5), log(2).
The null is conditional independence of Y and X given z/batch. The causal log
mean-ratio beta is known; the package's returned fold-change statistic is not
assumed to estimate this ratio without dilution. Eight untested background
Poisson(5) genes stabilize full-library QC summaries. Explicit `~z + batch`
avoids outcome-derived library covariates in the design. Standard QC nonzero
thresholds are 7; cellwise quantile ranges [0,1] retain all cells except assignment
eligibility. Every test is two-sided with the unchanged public family/resampling
policy. This evaluates small-family nominal testing, not genome-wide FDR.

The selected cases span low NT/complement permutations, high complement CRT
(empirical/skew-normal), high permutations/skew-normal, adjusted confounded CRT,
and deliberately omitted-z CRT. The last case is a **negative robustness
sentinel**, not an expected-green calibration case. Sparse Poisson(.02) responses
are all nulls and exercise QC rejection; tied constant-2 responses are all nulls
and use `~1` to stress discrete two-sided tails (including inherited issue237).
These two stress cases must not be pooled into the modelled-case calibration.
All-tied responses are a deliberately degenerate diagnostic, not NB data.

## Precision and interpretation

The independent unit is a regenerated dataset within a case. Endpoints sharing
cells are **never pooled as independent trials**. Summaries stratify by case,
response and nominal alpha (.01, .05, .1). Each summary reports valid tests,
QC failures, errors, invalid p-values, rejection fraction, Monte Carlo SE and
pointwise exact Clopper-Pearson 95% bounds. Invalid or failed tests are never
clipped/replaced. Conditional-on-valid rates and per-attempt rates are distinct.
The seed schedule is deterministic and stable under selecting/reordering cases;
DGP and test seeds are separated, and withr restores caller RNG state.

At true .05, 40 reps have MC SE .0345, 200 reps .0154, 1000 reps .0069 and
5000 reps .0031; a pilot cannot certify a nominal .05 rate within one percentage
point. Use 40 for bounded diagnosis, 200-1000 for follow-up, and 5000 per selected
case for precision near +/- .006 at .05. Approximate 95% halfwidth .01 needs about
1825 independent datasets at .05; use exact intervals in final interpretation.
Intervals are pointwise; multiple cases/tails are not simultaneous confidence
claims. Inspect invalid p-values/errors before interpreting rejection rates.
Do not threshold the whole suite into a green statistical-validity badge.

The manifest preserves git revision/dirty status, installed-file and source-file SHA256 hashes,
harness/setup hashes, settings, RNG kind/seed rule, R session and terminal time.
Each replicate preserves input SHA256, true effects, treated count, p/QC/stage,
warnings/errors and elapsed time. CSV checkpoints survive interruptions; only
`finished_utc` in manifest confirms terminal completion. Retain exact installed
source/archive identity alongside receipts; the checkout revision alone cannot
certify the binary used. No external data is needed or changed.

---

# Experimental Wilcoxon scientific pilot

From the current repository root with the exact candidate package loaded:

```sh
Rscript tools/examples/wilcoxon.R
Rscript tools/validation/wilcoxon-pilot.R 200 /tmp/wilcoxon-pilot
```

This produces actual rank-test p-values for independent 400-cell multinomial
screens, each with 1000 UMIs per cell, two targeting guides and two NT guides.
The first response has probability .01 in both arms under the exchangeable null;
the second has .02 in both arms, or .04 in treatment under the effect scenario.
Changes are balanced by a reservoir response, preserving total library size.
Replicates use seed `100000 + replicate`, paired across scenarios.

The confounded scenario changes the first response's probability to .04 in
batch 1 while treatment has 90% batch-1 cells and controls 10%. It has no causal
treatment effect on that response conditional on batch, but violates the marginal
rank-test null. Its rejections illustrate the lack of covariate adjustment;
they are not type-I error measurements under the rank-test null. The formula
includes batch to demonstrate that this alternative ignores it.

`summary.csv` reports raw p < .05 rates with 95% Wilson intervals over independent
screens. With 200 screens, uncertainty at .05 is substantial (standard error
about .015). `pvalues.csv` also retains two-pair BH .1 calls; raw rejection rates
are not family FDR estimates. `settings.rds` stores code hashes and session
information. Keep these identities with results. Do not present this pilot as
real-screen, donor-level, small-sample, high-MOI, backed-storage or broad FDR
validation. It does not establish equivalence to covariate-adjusted SCEPTRE.

`tests/testthat/test-wilcoxon.R` separately checks numerical parity with base R,
QC/index alignment, union/singleton/Bonferroni integration, complete-family
correction, all ties, missing-QC results, provenance, method transitions,
unsupported-mode guards, ownership, RNG preservation and save/reopen.
Reference agreement is an implementation check, not a calibration study.

See [the public API help](../../man/run_wilcoxon.Rd) and
[R's Wilcoxon documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/wilcox.test.html).
Distributed/Nextflow, backed responses and parallel execution are deliberately
unsupported. Follow-up requires method-aware execution manifests, matching
worker/collector provenance, direct/distributed tests and the same scientific
qualification before expanding those modes.

## Receipt boundary

Use a new output directory outside the clean checkout; existing paths are refused
before writes. `status.rds` starts as running, checkpoints progress, and becomes
complete only after all600scenario-runs (200independent screens per scenario, with paired seeds) and artifact hashes exist.
Errors/interrupts produce failed/interrupted status when R can handle them; a
hard kill leaves a running/stale checkpoint, never a success receipt. Do not
interpret an output directory or partial CSV as a completed run.

`settings.rds` binds clean Git HEAD and source-file hashes to the installed
package directory plus its actual loaded R/native build fingerprint. During
archive qualification set `SCEPTRE_PILOT_ARCHIVE` to that exact source archive;
its SHA256 is recorded too. The runtime and source identities are both retained,
not inferred from package version alone. Every screen must have both expected
responses, pass QC, and produce finite nonmissing probabilities in[0,1]; an
unexpected eligibility or probability result stops the run with failed status.
