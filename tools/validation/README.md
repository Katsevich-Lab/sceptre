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
600), comma-separated case IDs (default all), optional candidate setup R file.
Existing output directories are rejected. A setup file can set the option
`sceptre.validation.configure` to a function of `(object, case)` returning an
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
