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

The final validation receipt will record the immutable source commit, ordinary
installed/package checks with the oracle unset, five current-source probes,
manual old/new comparison and hosted CI. Preliminary failed test-authoring and
render attempts remain in local logs; they did not change production behavior.
