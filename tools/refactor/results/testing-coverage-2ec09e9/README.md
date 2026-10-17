# Testing before and after the refactor

This is a one-off developer measurement, not a CI step or package dependency.
Original 21f9ea0 is measured with its original package tests; 2ec09e9 is measured
with the permanent current package tests after the CI/manual-suite separation.
Both use the same R 4.6.0/macOS arm64/dependency environment and guarded ondisc.

- Before: 3,336 / 4,722 executable R lines, 70.6480%; 2,578 assertions / 56 blocks.
- After: 4,097 / 4,868 executable R lines, 84.1619%; 4,199 assertions / 138 blocks.
- Change: +13.5138 percentage points. Both test runs pass with no warnings/skips.

Coverage uses covr 3.6.5 package_coverage(type="tests") and tally_coverage(by="line").
Empty covr.flags disables native instrumentation; logs confirm ordinary native
compilation and no gcno files. Generated Rcpp R wrappers use covr's standard
exclusion. No C++ percentage, branch coverage or whole-workflow correctness
claim follows. Examples, vignettes, external native probes, old/new comparisons
and downstream tests do not contribute to this percentage.

The aggregate summary, per-file and per-line CSVs are readable directly. The raw
archive retains the exact measurement runner, full source/test SHA256 manifests,
coverage objects and install/test logs. Each archive member was read back and
hash-verified. It omits duplicate Git source archives and installed libraries.
Recreate source archives from the recorded commits; the retained runner includes
the private dependency-library location used for this measurement.

The assertion totals include 2,225 mock-data-helper assertions on each side;
counts alone do not describe production behavior coverage. Native mathematics,
controlled failure probes, manual consistency and actual consumer qualification
remain separately reported. The measured R gap includes the generic object
display/print/plot methods in s4_helpers.R and some validation branches.
