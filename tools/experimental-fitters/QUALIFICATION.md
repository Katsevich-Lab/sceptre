# Accelerated logistic fitter: bounded experimental qualification

The opt-in `fast_logistic` fitter is implemented and locally checked. `glm.fit`
remains the reference default. This is the same unpenalized likelihood with
numerically different fits and potentially different Monte Carlo draws/p-values.
It is **not** a bit-identical replacement or general calibration certification.
SPA, GPU, Python parity and other response models are not enabled.

## Exact artifacts and terminal checks

- Original component-model base: `840f33436a6811a45d712b78bc587d8a0bb8a25b`.
- Final tested source: `00a04f2e4b05f6f5a033e6f411e563fca89e83b8`, including
  separately admitted W02 empirical-cap and observed-score-guard foundations.
- Isolated Mac arm64 R 4.6.0 library: `/private/tmp/w07-library-repaired`.
  RcppEigen header version: **0.3.4.0.2**. New build dependency; old images do not
  qualify this code.
- Tested source archive SHA256:
  `cb8e9e55b5a3f8b6f279cfa34d999ae2e4b156986204a81e21adc5a5acf430f9`.
- Native DLL SHA256:
  `ae9ff3ccb3eaee87dfb12780d06443d7c0e7bb60ba429f1fb5950b2f4e970a05`.
- Targeted fitter tests: **77 assertions**, zero failures/errors/warnings/skips.
- Full package suite: **266 cases / 5,937 assertions**, zero
  failures/errors/warnings/skips. Runnable public help example passed.
- Physical old schema-three migration preserves current result identities and
  caller bytes. Removing only the new empty selector leaves the old spec exactly
  equal. Default snapshots/certificates remain compatible; candidate caches have
  their own method/native identity. Cold/warm results and RNG state agree.
- Publication pipeline adapter `41525a088b72920b3b30d554559f5820b01b95cf` passed
  actual candidate create/read/reference-stage rejection, genuine historical
  schema-three normalization, caller/input immutability and runtime-lock
  rejection. Adapter/test hashes matched before/after execution. This runtime
  lacked Wilcoxon APIs; it does not qualify their adapter lane.
- BiocCheck **1.49.30**, Bioconductor **3.24**, original submission wrapper flags
  `no-check-R-ver` and `no-check-bioc-help`: **0 errors / 1 inherited warning /
  3 inherited notes**. Strict category, individual function/line and display-call
  gate passed against the actual submitted source: **zero new instances**.
  This is the matched-submission lane, not an all-default-environment check.

Receipts are local at `/private/tmp/w07-repaired-gates`. Prior attempts remain at
`/private/tmp/w07-qualified-run*`: one missing-analysis test-harness error, then
one intentional additive-layout assertion failure and actual checker notes were
repaired, not hidden. The final full suite and gate above are fresh terminal
results on the repaired installed source.

## Independent-count pilot, with precision limits

W03's generator/truth was reused unchanged: 600 cells, empirical CRT,
skew-tail CRT and confounded CRT; 40 independent datasets per case; paired
reference/candidate methods. There were 240 actual public analyses / 720 endpoint
rows, with no errors, QC failures or invalid p-values. Independent datasets,
not endpoint rows, are the replication unit.

At raw alpha .05, both configurations rejected the known null **1/40** times
in empirical and skew cases (exact 95% interval **0.000633–0.131586**) and
**3/40** in the confounded case (**0.015742–0.203865**). Up/down endpoints were
detected **40/40** per case (exact lower 95% bound **0.911903**). These intervals
cannot rule out consequential miscalibration. They are pilot evidence, not
an FDR guarantee or a reason to replace the reference default.

Ten additional independent datasets each ran 12 dependent public calibration
pairs under each configuration. Dataset-level rejection fractions matched.
Those pairs must not be pooled as independent binomial trials. Calibration
provenance identifies the requested Newton/reference-fallback hybrid;
per-group transient fallback diagnostics are not persisted.

The pilot's exact implementation was `5cf2155907644f26da8feb5335d52336942f3b1c`
(test/example fix over the release scientific source). On repaired `00a04f2`,
all **13 tracked scientific function-body hashes** and the **native DLL hash**
match that measurement runtime. Only error-condition construction, help and
migration assertions changed. This does not relabel the original pilot as a
new integrated-package run.

### Why p-values can differ

All 360 paired endpoint input hashes match. Three p-values differ, all in
skew-CRT replicate 28: null **0.916 vs 0.888**, down **4.206499e-9 vs
2.512401e-9**, up **9.104598e-23 vs 3.137349e-21**. Recorded stages and rejection
counts agree; no general decision invariance is claimed.

A focused actual-bank diagnostic on repaired `00a04f2` found:

- Both fits converged without warnings: Newton 3 iterations, reference 4.
- Maximum fitted-probability difference **2.15137e-9**, coefficient difference
  **9.47666e-9**. The reference equals the tight glm fit here.
- Under the same 5,498 draws and fixed native seed, materialized draw-bank
  hashes differ and all 5,498 draw rows differ. Observed scores and shared
  response-score inputs are identical; null-score hashes differ and low-level
  results reproduce the three public p-value differences.

This establishes a fitted-probability-to-sampled-bank difference in this
fixture. The precise internal Boost branch causing stream divergence was not
localized. No tolerance was tuned to reproduce reference p-values.

## Runtime, allocation and real-data scope

Twenty fit-only timing repetitions, six design columns:

| Cells | Reference cold | Candidate cold | Candidate prepared |
|---:|---:|---:|---:|
| 5,000 | 4 ms | 1 ms | 1 ms |
| 20,000 | 15 ms | 4 ms | 3 ms |

Warm reconstruction was below the millisecond timer resolution, not zero cost.
Tracked R allocation volume fell from **9,154,464 to 759,608 bytes** and from
**36,574,464 to 3,028,664 bytes** respectively. `Rprofmem` excludes Eigen/native
malloc: these are **not total-memory ratios**. Entire pilot peak RSS was
422,215,680 bytes including libraries and fixtures, not a fitter-specific peak.
Compact fits were 2,440 candidate versus 1,904 reference bytes. Tight-glm
probability errors were at most 2.82e-13 in these two fit fixtures. No
end-to-end, high-dimensional, full-screen or cluster speed claim follows.

A read-only targeted-panel DC-TAP comparison used 3,000 evenly spaced cells,
six illustrative threshold-5 guide groups and two responses. Original observed
panel and guide totals were computed before cell subsetting. All 12 seeded
p-values/fold changes agree; maximum probability difference was 2.34121e-9 and
coefficient difference 6.97668e-8. This supplies numerical comparisons, not real
null/effect truth, whole-transcriptome depth or biological discovery claims.
All real-derived tables, cell selections and manifests stay **local**, outside
published Git ancestry. Only these compact aggregate metrics are recorded here.

## Reproduction and remaining limits

From an isolated installation and this repository root, use:

```sh
R_LIBS=/isolated/library Rscript tools/experimental-fitters/example.R
R_LIBS=/isolated/library Rscript tools/experimental-fitters/evaluate.R \
  /new/evidence /path/to/tools/validation/scientific-harness.R
R_LIBS=/isolated/library Rscript tools/experimental-fitters/diagnose-replicate.R \
  /new/replicate28.rds /path/to/tools/validation/scientific-harness.R
```

Use the shared heavy lease for execution. The evaluator deliberately refuses an
existing evidence directory. Preserve source/installed/harness hashes and
source provenance. See `SOURCE.md` for pinned upstream source.

W03 independently approved restricted opt-in experimental method/sampler
admission after source mathematics, actual bank/output reconciliation and
terminal receipt readback. This closes that bounded review gate; W03 did not
rerun the full suite or checker. Larger statistical studies, broad estimability,
Eigen allocation/OOM/interrupt stress, floor/cross-platform and container
qualification, and exact final integrated-runtime checks remain separate gates. Current observed-score safety is a precise guard, not broad
estimability proof. Keep this candidate opt-in and direct-R only; do not use
historical SPA benchmark summaries as qualification for it.
