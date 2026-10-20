# Installed-package coverage

Run from a dedicated clean checkout, with the package's test dependencies plus
`covr`, `jsonlite` and `xml2` installed:

```sh
Rscript --vanilla tools/coverage/measure.R . /absolute/fresh/report-directory
```

The runner installs an instrumented package into `report-directory/library`,
then runs its **installed** test suite. It never installs into ordinary user
libraries. Builds use at most two compiler jobs; BLAS threads are limited to one.
Local overnight runs must first acquire the program's shared resource lease.
Do not edit source/tests or the runner while it executes. The runner binds input
hashes before compilation and rejects a changed source/test tree.

Open `index.html` in a browser. Each file expands to annotated source: red lines
were not reached, blue lines were reached. The artifact also contains per-line
and per-file CSVs, `summary.json`, Cobertura XML, raw covr RDS, session information,
and separate measured R/native SVG badges. No coverage-service account is needed.
The workflow's job summary links the artifact conceptually; GitHub requires
opening the run's Artifacts panel, downloading it, and opening the HTML locally.
No PR comment or public hosting is required, including for fork PRs.

## What the metric means

The summary counts unique physical executable lines reported by covr/gcov, using
the minimum hit count when several instrumented expressions share a line.
R and native counts are always separate. This is **not branch coverage**, a
scientific validation score, or proof of native memory safety. Tests run native
code through the installed shared library, but standalone sourceCpp bridges and
external native probes are not included in this measurement. Examples and
vignettes are not executed by this tests-only lane.

Only covr's standard generated wrappers are excluded (`R/RcppExports.R`,
`src/RcppExports.cpp`; no handwritten exclusions are added). Source inventory
files without executable records are explicitly listed. Declaration/documentation
files may legitimately have no executable records; an unavailable native toolchain
must never silently become 100% native coverage. CI requires at least one measured
native line and fails if instrumentation yields none.

GCC/gcov is used on Linux; Apple clang uses the matching `llvm-cov gcov` interface
on macOS. Compiler flags and tool paths are retained. This macOS baseline used
Apple clang21, R4.6.0 arm64, covr3.6.5, and locally installed ondisc1.2.0; it is not
proof of the pinned ondisc commit. CI explicitly installs/verifies the existing
qualified ondisc source pin. Unexecuted overflow/capacity guards remain visible;
forcing impossible allocations merely to raise this percentage is inappropriate.

## Review and badges

The new `coverage.yaml` workflow runs on PRs, main/dev pushes and manual dispatch,
with read-only contents permission and checkout credentials disabled. It uploads
annotated HTML plus raw evidence/logs even when later steps fail. Metrics badges
are generated from that run's measured values and belong to the same artifact;
they are not hardcoded promises about main or an untested integration revision.
A README workflow badge can link the actual workflow run status after W10 publishes
the branch. Before publication its absence/unknown state must remain visible.
There is no external coverage-service token configured or required. A permanent
public numeric badge URL would require separately authorized artifact hosting or
service setup; this work does not silently publish files or invent that service.

The pipeline needs its own measured CI/nf-test evidence. An R coverage number
cannot stand in for pipeline coverage. No pipeline files are changed by W01.

## Baseline and gap-driven tests

At package base `840f33436a6811a45d712b78bc587d8a0bb8a25b`, the installed suite
reported **5641 pass, 0 fail, 0 warn, 0 skip**. R coverage was **7046/7668 = 91.8884%**;
native was **953/970 = 98.2474%**. This is a local baseline, not hosted CI.
The first runner invocation returned1 after writing the complete coverage outputs:
its trailing source was edited while R was reading it, causing a parse error.
That tooling error is preserved separately from the passing installed suite.
The final frozen-run receipt supersedes it as verification of the runner itself.

Uncovered native lines were mainly capacity/overflow guards in sampling/count
code, a tail-fit monotonicity rejection, mixture component relabeling, and an
invalid combination-count guard. `inference.h` contains declarations only;
`sampling_safety.h`, `result_lifetime.h`, and `pair_eligibility.h` have measured
records. Raw native gcov output is retained with the local baseline evidence.

The lowest R coverage included getters/setters, validation, backed IO and plotting.
Material lifecycle gaps led to new tests for guide-count replacement: independently
recomputed full-library count/maximum-guide summaries, invalidation of downstream
preparation/fits, unchanged response counts, input ownership, and rejection of
misaligned cell/feature labels. Cache tests cover duplicate/unnamed update
rejection and explicit invalidation/removal without mutating the old cache.

Additional hand-calculated result checks cover QC-aware guide Bonferroni counts,
all-failed-QC summaries, the p=1 cap, calibration retaining guide rows, finite-family
BH/Bonferroni correction, scalar alpha despite a same-named tail column, and a
strict binary-exact significance boundary. Two deliberately wrong significance
vectors are proved to fail the same assertion used on genuine output. Existing
`tools/refactor/check-test-strength.R` also checks corrupted workflow results and
wrong/constant/correlated/repeated native draws; it is retained unchanged.

These are deterministic regression/algebra checks, not operating-characteristic
calibration. W03 owns independent simulation validation; W02 owns native safety.
