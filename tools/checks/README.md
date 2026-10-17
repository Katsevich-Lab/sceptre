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
