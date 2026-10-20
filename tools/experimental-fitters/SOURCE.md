# Source and acceptance boundary

The native logistic solver and R prepared-design helpers are reused from
ZiangNiu6/sceptre, `codex/tail-fitting-study`, verified by authenticated `gh` at
`7ea1aefd246f87965a16a1cfc5a4bf5b1901e05a` on 2026-10-19 US Eastern.

- `src/grna_logistic_fit.cpp`: original numerical solver, GPL-3; integration adds
  attribution and rejects malformed zero-column/empty prepared payloads before
  intercept access and widens column-offset products before multiplication.
  No numerical algorithm rewrite.
- `R/grna_logistic_fit.R`: prepared design and index helpers from the same source.
- Integration, selector, certificates and tests are new for the schema-3 base.

Source URL: https://github.com/ZiangNiu6/sceptre/tree/7ea1aefd246f87965a16a1cfc5a4bf5b1901e05a

This is a same-likelihood, numerically different fitter. It is not a bit-identical
replacement for glm.fit. Native allocation/interrupt stress for the new Eigen
solver has not been established by the existing reference-kernel safety suite. The solver conditions a leading-intercept design and
uses Newton updates with monotone likelihood checks and score/deviance stopping.
Failures return to reference glm.fit; they do not select a different scientific
model. Candidate probabilities are reconstructed through the existing reference
inverse link so cold and compact-cache warm use agree exactly.

Build dependency: RcppEigen headers are required even when the selector stays at
the default. Do not qualify old containers or archives for this code. The
pipeline config does not expose this selector. The package execution manifest
guard rejects candidate objects. The reference pipeline direct adapter also
requires its own explicit selector guard; that companion change is owned by W10
and is a prerequisite for integrating this feature. Until that guard is
qualified, do not feed candidate objects to a reference pipeline.

SPA is separate. The same upstream fork includes exact/fast CRT SPA and empirical
studentization variants. Live Compass's 2026-08-30 Gasperini summary reports up to
threefold p-value-stage and 30% total speedups for its CRT benchmark; it is a
derivative historical summary and does not establish validity for this revision,
permutation tests or this accelerated logistic fitter. No SPA source/routing is
integrated here. Independent SPA calibration, extreme-tail/root failure and
fallback qualification remain deferred.
