# Manual refactor consistency suite

This suite is deliberately outside package test discovery and excluded from
`R CMD build`. No GitHub Actions workflow invokes it, including scheduled or
manually triggered workflows. It compares the refactored package to original
commit `21f9ea098b69c444a884a49bc254de711968e3b5` on frozen inputs.

From a clean committed checkout, run:

```sh
Rscript --vanilla tools/refactor/compare.R /absolute/fresh/output
```

An optional second argument selects a committed candidate Git ref. The default
is HEAD; uncommitted edits are never included. Output must be outside the source
checkout. Run with the same dependency library you intend to validate, e.g. set
`R_LIBS` before launching if using a private ondisc installation. Developer tools
`pkgbuild`, `pkgload`, `testthat`, `digest` and `jsonlite` must already be available;
backed comparisons run when ondisc is installed.

The command archives both Git refs, performs a fresh release installation of the
candidate, runs verified original source in a separate process, and tests the
installed candidate against those observations in another process. Baseline
source/input hashes, original numerical environment and complete sampler banks
retain their existing checks. No candidate snapshots replace expected output.
Reports include `comparison.json`, per-contract `checks.csv`, detailed RDS and
full logs. Failed stages return nonzero and retain their logs. A predeclared
platform-specific archive-stream check can skip on another toolchain; the
same-environment source-stream comparison still runs.

`testthat/` contains the exact pre-separation helper/test bytes and complete
historical fixtures. `relocation.json` binds every preserved file to finishing
commit 04df2f0 and its SHA256. These files are not used by permanent CI. Only
shared mathematical inputs are also retained in package tests. Permanent tests
check mathematics, reproducibility, current sparse/backed and serial/fork parity,
cache invalidation and current saved-object roundtrips without this suite.

The other explicit refactor probes and actual pipeline/resource tools remain
available; this command is not their combined qualification. Previous evidence
and its measured source pins remain unchanged. After maintainer review, retain
a tagged source/evidence snapshot before deleting this subtree from active work.
Do not remove permanent tests or the lasting root architecture guide with it.
