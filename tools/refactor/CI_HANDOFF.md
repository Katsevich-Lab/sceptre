# CI and refactor review tooling separation

This follow-up changes tests, workflow wiring and developer documentation only.
Production `R/`, `src/`, `NAMESPACE` and `DESCRIPTION` remain byte-identical to
qualified source `0a9f9d879db6a52e3daef237daa1711b9d542976`.
The original acceptance/resource/consumer reports retain their measured pins.

| Material | Active home | Routine CI? | R source package? |
|---|---|---|---|
| Current behavior, indexing, mode parity, input nonmutation, cache and RDS contracts | `tests/testthat/` | Yes | Yes |
| Current mathematical and native ownership/exception probes | `tools/checks/` | Yes | No |
| Original package/output/sampler comparisons and complete old fixtures | `tools/refactor/consistency/` | No | No |
| One-command immutable old/new comparison | `tools/refactor/compare.R` | No | No |
| Tim's reading guide, historical consumer/resource runners and evidence | `tools/refactor/` | No | No |
| Lasting module and coordinate guide | root `ARCHITECTURE.md` | Documentation | No |

The permanent suite no longer fetches/builds old SCEPTRE, reads the historical
API/output snapshots, compiles a frozen old sampler, or requires a runtime
oracle. Mathematical input cases are retained. Aggregation expectations use
dense column sums and explicit legacy shape/indexing rules, rather than a
saved old implementation. Complete current workflows are checked for
reproducibility; current execution/storage modes are compared with one another.
Historical equality remains explicit in the manual suite, rather than being
quietly dropped or replaced by a disabled CI test.

The original helper, three comparison test files, full fixtures and aggregation
reference are preserved byte-for-byte from finishing commit 04df2f0. The 32-file
SHA256 mapping is in `consistency/relocation.json`; no historical expectation was
regenerated. Historical mixed native probes also remain available outside CI.

CI retains its Linux/macOS/Windows package matrix, Air and ondisc-free lane.
The exact guarded ondisc dependency pin remains necessary for backed safety;
it is separate from old-SCEPTRE comparisons. Both `tools/checks/` and
`tools/refactor/` are excluded by `.Rbuildignore`. Branch-specific duplicate push
runs are removed; the refactor branch is checked through its pull request.

For a reproducible comparison, see `consistency/README.md`. For a guided source
review, read `walkthrough.qmd`. Tim has not been messaged or assigned work by
this change. The tooling stays available until maintainer review; later deletion
should follow a tagged source/evidence snapshot.

## Validation record

Local validation is retained in [the evidence receipt](results/ci-handoff-370897a/qualification.json)
and [restorable raw logs/results](results/ci-handoff-370897a/README.md).

| Check | Qualified source | Outcome |
|---|---|---|
| Installed permanent suite; `CI=true`, oracle unset, guarded ondisc required | 708f37a | 4,199/4,199 expectations; 138 blocks; zero failures/errors/warnings/skips |
| Current native mathematics | 708f37a | 9 independent expectations |
| Controlled current sampler/count/dispersion/result probes | 708f37a | 9 / 214 / 14 / 959 cases; negative controls retained |
| Full `R CMD build` and `R CMD check --no-manual` | 708f37a | Status OK; examples and rebuilt vignettes; no historical oracle |
| Source tarball inspection | 708f37a | Both developer directories and old output/source oracles absent |
| Front-door old/new comparison | 370897a vs original 21f9ea0 | 1,265/1,265 expectations; 31 blocks; zero failures/errors/warnings/skips |
| Historical relocation | 04df2f0 to new paths | All 32 files byte-exact; SHA256 verified |
| Initial walkthrough (before illustrated revision) | QMD hash in receipt | Rendered HTML; executed indexing result and live source snippets checked |

The receipt verifies identical permanent source/test/workflow bytes between
708f37a and 370897a. The latter changes only the manual wrapper and Git
attributes. Subsequent handoff evidence/documentation is build-excluded.
Production files remain identical to 0a9f9d8. No resource or actual-consumer
rerun is claimed for this test/tooling separation.

Hosted checks belong to the published PR head; consult
[PR 234 checks](https://github.com/Katsevich-Lab/sceptre/pull/234/checks).
The local receipt does not substitute for their outcome or claim future heads
have passed. Preliminary failed test-authoring, wrapper and render attempts
are retained separately in the evidence archive.

## Hosted result and guide revision

All six jobs pass for published code/test/workflow head `e230067` in
[run 38987592424](https://github.com/Katsevich-Lab/sceptre/actions/runs/38987592424).
The subsequent illustrated tour and bug appendix change developer documentation
and figures only. Their renders and local links are checked separately; package
source, tests and CI inputs retain the receipt's exact bytes.

The guide revision's source/HTML hashes, executed example, links and diagram
checks are recorded separately in [TOUR_VALIDATION.json](TOUR_VALIDATION.json).
