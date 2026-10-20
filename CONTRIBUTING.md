# Contributing to sceptre

Contributions to `sceptre` are welcome, including code, tests, and documentation.

Participation in the project is subject to the [Code of Conduct](CODE_OF_CONDUCT.md).
Please do not include private data, credentials, protected health information,
or other sensitive information in pull requests, examples, or test fixtures.

For bug reports, feature requests, and questions about using `sceptre`, see
[Getting help with sceptre](SUPPORT.md).

## Contributing code and documentation

Small bug fixes and documentation corrections may be submitted as pull requests
without prior discussion.

Before implementing a substantial feature, public API change, or change to
statistical methodology, start a [feature-request discussion](https://github.com/Katsevich-Lab/sceptre/discussions/categories/feature-requests).
Early discussion helps the maintainers assess scope, compatibility, and
validation requirements before substantial work is invested. Maintainers may
promote a sufficiently scoped proposal to an issue for implementation tracking.

### Code, tests, and documentation

For the code's main parts, cell numbering, numerical calculations and pipeline
compatibility, see [Internal architecture](ARCHITECTURE.md).

| Material | Where it belongs | How long it stays |
|---|---|---|
| Package implementation, help, examples and datasets | `R/`, `src/`, `man/`, `vignettes/`, `data/`, `inst/` | Permanent; maintain with the package |
| Tests of current behavior and bug fixes | `tests/testthat/` | Permanent, including regression tests first added during the refactor |
| CI and current-source C++ checks | `.github/workflows/`, [tools/checks/](tools/checks/README.md), `air.toml` | Permanent maintainer tooling |
| Test-input recipes and provenance | `tools/fixtures/`, `tests/testthat/fixtures/README.md` | Permanent; never regenerate expected scientific outputs in CI |
| Module-map drawing and editable source | `tools/architecture/` | Permanent; maintain with `ARCHITECTURE.md` |
| Current architecture and contributor/project guidance | Root `ARCHITECTURE.md`, this guide, `README`, `ROADMAP.md`, `GOVERNANCE.md` and other project policies | Permanent; keep current |
| Refactor tour, manual comparison and validation receipt | [tools/refactor/](tools/refactor/README.md); its published archive retains the specification, plans and bulk evidence | Temporary; retain for maintainer review, then archive before retiring |

**Build exclusion does not mean temporary.** `.Rbuildignore` keeps maintainer
material out of the distributed R source package. Permanent `tools/checks/`,
`tools/fixtures/` and `tools/architecture/`, and temporary `tools/refactor/`, are
excluded, as are several lasting root guides.

Routine CI uses current-package tests and current-source checks. It must not
depend on `tools/refactor/`, an old SCEPTRE installation or historical comparison
outputs. Some permanent tests still have `refactor` in their filenames; that
records their origin, not a plan to delete them. Remove temporary review links
when retiring the refactor directory, and keep the lasting architecture guide
independent of it.

- Keep changes focused on a single problem or feature.
- Add or update tests for behavior changes and bug fixes.
- Format handwritten R code with Air as described below.
- Update user-facing documentation when public behavior changes.

Changes affecting model fitting, calibration, power, numerical results, or other
statistical behavior generally require more than unit tests. Describe the relevant
simulation, benchmark, or reproducible analysis in the pull request.

### Fresh development setup

Use a supported R version (the package declares R >= 4.1), a C++ toolchain,
and the dependencies in `DESCRIPTION`. Work on a focused branch from the
maintainers' requested base. To keep development dependencies and compiled
artifacts separate from another project's installation, use an isolated library:

```r
lib <- file.path(tempdir(), "sceptre-dev-library")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
# Install remotes into lib first if it is not available.
remotes::install_deps(".", dependencies = TRUE, lib = lib, upgrade = "never")
pkgload::load_all(".")
testthat::test_local(".", filter = "wilcoxon")
```

`pkgload::load_all()` compiles current native sources; an installed package with
an older DLL is not a substitute for checking a changed native implementation.
For backed-data tests, follow the exact `ondisc` pin and system dependencies in
`.github/workflows/check-standard.yaml`. Keep it optional: absence must not break
in-memory use. Never regenerate scientific fixtures to make a failing test pass.

### Choose a validation tier

Run the smallest meaningful checks while developing, then the applicable wider
checks before review. Record commands, source revision, platform, dependency
versions, skipped lanes, and terminal results; a submitted job is not a pass.

| Tier | Command or entry point | What it establishes |
|---|---|---|
| Fast | `testthat::test_local(".", filter = "wilcoxon")` (replace the filter for your change) | Focused behavior and regression checks using current source |
| Full package | `testthat::test_local(".")`, then `devtools::check()` | Package-wide tests, examples, documentation and installation checks; inspect skips |
| Native | `Rscript --vanilla tools/checks/test-native-components.R .` and the ownership/exception probes listed in `tools/checks/README.md` | Mathematical components and declared error/lifetime boundaries; these compile temporary instrumentation |
| Scientific | `Rscript tools/validation/wilcoxon-pilot.R 200 /tmp/wilcoxon-pilot` plus the applicable method-specific harness | Actual simulated null/effect p-values with Monte Carlo uncertainty; reference parity alone is not calibration |
| Cross-surface | Same exact archive/runtime and complete family in direct and distributed workflows | Supported-mode equality and companion pipeline compatibility; an unsupported method must reject before execution |

Statistical changes need an independent oracle or known-data-generating process,
not just snapshots produced by the implementation. Show known nulls and effects,
confounding and degenerate cases, conditioning and replication units, and limits
of precision. Test that the regression would fail for the incorrect result.
Keep larger scientific runs reproducible outside routine tests. Bound local CPU,
memory and runtime, and coordinate shared resources before starting heavy work.

For Bioconductor qualification, compare both diagnostic categories and flagged
instances against the actual submission with the same checker/options. Do not
suppress new diagnostics, weaken checks, or treat a different checker environment
as a matched baseline. Retain environmental diagnostics separately.

### Worked extension: an alternative test or fitter

`R/wilcoxon.R` is a concrete experimental test extension, with its public help in
`man/run_wilcoxon.Rd`, regressions in `tests/testthat/test-wilcoxon.R`, and a
runnable example in `tools/examples/wilcoxon.R`. Read it alongside the execution
and lifecycle sections of `ARCHITECTURE.md` before adding another method.

1. Define the null, estimand, response scale, control population, independence
   assumptions and supported modes first. Wilcoxon explicitly uses unadjusted
   normalized ranks, not SCEPTRE's covariate-adjusted test. A descriptive mean
   ratio is labelled separately from the rank estimand and has no effect SE.
2. Reuse prepared cells, `CellView`, QC and storage access. Retain full-library
   depth covariates; do not recompute them from requested genes. Preserve failed
   QC rows and multi-guide semantics.
3. Assemble all requested pairs and apply guide integration and correction once
   over the complete family. A worker's subset is not an independent family.
   Commit only after plan/family/probability validation succeeds.
4. Save actual method, scale, fitting policy and execution provenance. Keep method
   transitions and caches honest. An accelerated fitter must bind its method to
   cache identity and specify convergence/failure behavior. A different fitter
   is not automatically a bit-identical implementation of the old one.
5. Add independent numerical checks, actual null/effect simulations, failure and
   unsupported-mode tests, save/reopen and nonmutation checks. Preserve default
   RNG/draw ordering. Use explicit opt-in until scientific acceptance is complete.
6. Update help, examples, reports and companion pipeline/runtime contracts.
   `run_wilcoxon()` currently supports direct serial in-memory low-MOI NT controls only;
   distributed/Nextflow support is deferred and guarded. Do not advertise parity
   for that lane. See its example for method comparison using separate objects.

The experimental Wilcoxon path delegates numerical p-values to the installed R
stats version. Numerical cancellation or underflow can produce zero even when
the tail probability is representable. Raw zeros remain in results; log-scale
plots use a display floor only, without repairing probabilities. Extreme-tail
magnitudes and extreme-tail calibration remain unqualified. Reference agreement
with base R does not independently validate these tails.

Exported functions documented in `man/` are the public interface. Component
slots, `RunPlan`, views, fitters, caches and `sceptre:::` helpers are private;
examples for contributors do not make these stable user APIs. The pipeline has
specific internal contracts documented in `ARCHITECTURE.md` and
`tools/hardening/DISTRIBUTED_EXECUTION_CONTRACT.md`; coordinate changes with its
maintainers rather than adding incidental private calls. Do not use compatibility
aliases as new storage. Public API and methodological changes still require the
maintainer discussion described above.

Regenerate help with `roxygen2::roxygenise()` using the repository's roxygen
version, inspect the generated diff, and build/check examples. Avoid unrelated
regeneration of help or bindings. Use Air for handwritten R as below; native code
and test fixtures need their own review.

### Formatting R code

Use Air 0.9.0 as the single R formatter. The repository's `air.toml` specifies
four-space indentation and a preferred width of 80 characters. Format only the
handwritten package source, then check it with the same command used by CI:

```sh
air format R
air format --check R
```

Keep Air's default exclusions for generated bindings. Immutable fixtures and
`tools/` (including reference code and retained results) are explicitly excluded;
do not recursively format tests or reference/evidence directories.

Air does not guarantee that every line is at most 80 characters. Wrap ordinary
prose comments and documentation separately where useful, preserving code
examples, fenced blocks, URLs, mathematical notation and executable strings.
BiocCheck remains a separate Bioconductor checker; Air does not replace its
checks or guarantee that every style note is resolved.

### Before submitting a pull request

Before submitting a pull request, run a package check:

```r
devtools::check()
```

Open the pull request against `main` and explain:

- the problem being addressed;
- the proposed change;
- the tests and checks performed;
- relevant documentation changes; and
- any statistical, compatibility, or performance implications.

Also verify in the review description:

- Public changes have help, an executable example, supported-mode limits and
  explicit method provenance.
- Tests use current source and meaningfully detect the wrong result; skips and
  precision limits are stated.
- Failed QC, complete-family correction, object ownership and method/cache
  transitions remain correct.
- Defaults, RNG, native behavior and scientific fixtures are preserved unless an
  intentional numerical correction is separately justified.
- Applicable native, scientific, formatter, package and matched submission checks
  have terminal receipts; pipeline companion changes use the same exact runtime.
- No private data or unrelated files/libraries are included.

A maintainer will review the contribution and may request revisions, additional
tests, or further validation before merging. Maintainers may decline changes
that fall outside the project's scope or current roadmap.

## License

Contributions accepted into the repository are distributed under the project's
[GPL-3 license](LICENSE.md).
