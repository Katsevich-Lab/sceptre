# Independent calibration/power pilot, October 19, 2026

The pilot completed **360 independent simulated cell datasets** (nine cases,
40 per case, 600 cells per dataset). It found a reproducible inherited failure:
**40/40 constant-count datasets abort at the completed-result probability guard**.
The six modelled NB regimes returned no errors, QC exclusions or invalid p-values;
their observed null rejection fractions do not provide precise calibration
certification at this replication count. No default scientific code changed.

## Exact execution and identity

- Pinned package R/native input: `840f33436a6811a45d712b78bc587d8a0bb8a25b`.
- Source checkout/harness at execution: `d53166030b7784447a768d64d3bee385e046671e`.
  Changes since the pinned package are excluded tooling plus one actual-count test.
- Isolated installation: `/private/tmp/w03-scientific-library/sceptre`.
  `install.log` records successful source compilation/installation. No ordinary
  user library was written. `manifest.rds` contains the per-file installed-binary
  and source R/native SHA256 identities and exact R/dependency session.
- Command, under the shared atomic local-heavy lease, with BLAS/OpenMP=1,
  MAKEFLAGS=-j2 and data.table threads=1:

  ```sh
  R_LIBS=/private/tmp/w03-scientific-library Rscript tools/validation/run-scientific.R \
    /private/tmp/w03-scientific-pilot 40 600 all
  ```

- UTC start `2026-10-20 01:40:18`, terminal `2026-10-20 01:40:53` (October19 local).
- `replicates.csv` SHA256:
  `c8fdbdee95e0d05fb7bc2e8a502981504daf26bb15110cfe7f236bae9f6d0d89`.
- No setup/candidate method was used. No source changed during the running job.
  Subsequent post-QC-hook and namespace qualification edits are not relabelled
  as part of this run; they do not change default DGP/test computation.

## Results at alpha .05

Intervals are **pointwise exact binomial 95% intervals**, separately for each
response/case. Outcomes sharing cells are never pooled as independent trials.
All modelled-case null endpoints have true beta=0; their two effect endpoints
have beta=log(.5) and beta=log(2). See README and the saved case manifest for the
complete DGP and unchanged public inference policy.

| Case | Null rejections / valid datasets | Rate | Exact 95% interval |
|---|---:|---:|---:|
| Low MOI NT empirical permutation | 3/40 | .075 | [.0157, .2039] |
| Low MOI complement empirical permutation | 3/40 | .075 | [.0157, .2039] |
| High MOI empirical CRT | 1/40 | .025 | [.0006, .1316] |
| High MOI skew-normal CRT | 1/40 | .025 | [.0006, .1316] |
| High MOI skew-normal permutation | 0/40 | 0 | [0, .0881] |
| Confounded, z-adjusted empirical CRT | 3/40 | .075 | [.0157, .2039] |
| Confounded CRT, z deliberately omitted | 40/40 | 1 | [.9119, 1] |

In each of the six modelled regimes, **both effects were detected in40/40**
datasets (separately, each exact95% interval [.9119,1]). The intentionally
misspecified case detected the down endpoint in4/40 and up in40/40, illustrating
how positive confounding can mask a negative effect as well as cause false
positives. This is not an expected-green case or a method-correctness failure
under its assumptions. The six adjusted/randomized cases are small pilots and
cannot establish equivalence between methods, precise type-I error, or FDR.

Sparse Poisson(.02) responses were all null: each endpoint failed nonzero QC
in40/40 datasets. They are eligibility exclusions, not p=1 observations and not
calibration successes.

For constant-2 responses under an intercept-only design, all40 datasets failed
with `Completed result p-values must be numeric probabilities or missing.`
The public workflow exposes **no returned raw invalid p-value**; all three
requested endpoint records contain that workflow error. W02's independent native
oracle identifies the inherited two-sided all-ties p=2. These are two distinct
pieces of evidence: actual-count workflow failure here, native tail diagnosis
there. The existing guard prevents saving an invalid result but does not make
the analysis complete successfully.

## Deliberate issue237 correction assessment

W03 independently inspected W02's actual `src/inference_tail.cpp` diff in
`e33ee20ad71debec92fa09bba5a1bf7046675aa0`: an outer `min(1,...)` around the existing
inclusive doubled smaller empirical tail, plus an explanatory comment. The fix
is scientifically justified as an explicit probability-domain correction.
One-sided tails, tail counts, loops, scores and random draws do not change.
Only values above1 change, so decisions at any alpha<1 and the .02 stage
threshold cannot change. BH decisions at alpha<1 also cannot change: those old
values cannot pass a BH threshold, capped1 values still cannot, and ranks below1
remain unchanged. W02 reports195 terminal oracle/RNG/stage/input checks; that
receipt is attributed to W02, not counted as W03 execution. W11 independently
reports5460 exhaustive small discrete native oracle checks, also separate.

**Admit this narrow correction**, then rerun `all_tied` against the exact final
integrated installation. This pilot intentionally retains the actual baseline
failure; it has not been relabelled as testing the corrected binary.

## Routine checks and limitations

`check-harness.R` passed generator/aggregation integrity checks. The installed
actual-count default-method regression recorded **4 passes,0 failures,0 warnings,
0 skips** (`sentinels.log`). This is a strong-effect direction/detection sentinel,
not a calibration estimate or comprehensive package test run.

At nominal .05,40 independent datasets have theoretical MC SE .0345; the displayed
exact bounds are the authority, including when plug-in MC SE is0. About1825
independent datasets per selected case are needed for approximate95% halfwidth
.01, and5000 gives MC SE .0031. The coordinator directed release of the lease
after this pilot; a larger batch follows the other workers' queued checks.
Do not combine repeated seeds from an expansion with the pilot as independent
replicates. The harness seed schedule deliberately allows exact replay.

Not assessed here: genome-wide FDR, assignment/mixture errors, multi-guide target
integration, clustered/donor dependence, rare subpopulation effects, real-screen
calibration, backed CRT, distributed runtime equality, alternative fitter or
Wilcoxon, and exact final integration binary. Existing reference-equality
receipts are not evidence for these claims. No BiocCheck or release qualification
was run by W03; W04/W10 own those gates. No live W03 job remains after the pilot.
