# Permanent source checks

These checks validate the current implementation without loading an old package,
old implementation, historical expected output or refactor evidence archive.
They run in routine CI alongside the self-contained package tests.

- `test-native-components.R`: direct effect/score/tail mathematics.
- `test-sampling-ownership.R`: controlled pointer-boundary failures.
- `test-count-ownership.R`: controlled count allocation/error paths and nonmutation.
- `test-dispersion-exceptions.R`: current numerical-versus-system exception policy.
- `test-result-ownership.R`: current result-buffer lifetime, controlled failures
  and negative control. Unfailed current calls are repeatability controls, not
  numerical-reference claims; package/component tests supply mathematical checks.

Each source probe compiles test-only instrumentation in a temporary directory.
No production failure-injection hook is added. These certify their declared
buffer/error boundaries, not every Rcpp allocation or actual operating-system OOM.
This directory is excluded from the source package by `.Rbuildignore`.

Historical equality checks, including the original forms of mixed probes, live
in `tools/refactor/` and have no CI entry point. See its `README.md`.

## Deliberate empirical p-value correction (issue #237)

The original two-sided empirical value doubled the smaller add-one, inclusive
tail without a bound. For null scores `c(-1, 0, 0, 2)` at observed zero, both
one-sided values are `4/5`, so it returned `8/5`; all ties returned `2`. The
corrected value is `min(1, 2 * min(left, right))`. This is a deliberate numerical
bug fix, separate from mechanical component refactoring. One-sided values and
the treatment of ties/add-one correction are unchanged. Empty draws now yield
one rather than two for the two-sided value.

Only original values above one change. Hence the stage threshold `.02` and all
individual rejection decisions at alpha below one stay unchanged. Capping also
preserves BH rejections at levels below one: original values above one cannot
pass any BH cutoff and stay above every such cutoff when capped. No sampler,
score, budget, selected draw segment, fit or random-number call changes.
The source probe checks independent tail counts and bounds, zero/one draws, both
extremes and mixed/all ties, nonmutation, unchanged R RNG and staged draw values.
Frozen historical equality checks that intentionally assert values above one
need this exception when comparing against pre-correction versions.
