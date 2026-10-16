# Refactor status — October 16, 2026

**Complete, with two explicitly accepted bounded memory costs.**
Production checkpoint: `0a9f9d879db6a52e3daef237daa1711b9d542976`.
Original authority: `21f9ea098b69c444a884a49bc254de711968e3b5`.
Branch: `codex/refactor`; [draft PR 234](https://github.com/Katsevich-Lab/sceptre/pull/234).
Nothing is merged or released.

The [final acceptance report](FINAL_ACCEPTANCE.md) binds source pins, measured
results, accepted costs and limits. The
[evidence directory](results/final-0a9f9d8/README.md) retains raw successful,
failed and invalidated attempts through verified content-addressed archives.

| Gate | Final outcome |
| --- | --- |
| Installed package at 0a | 4,740/4,740 checks, 154 blocks; zero failures/errors/warnings/skips; backed coverage required |
| Native at 0a | 112 component expectations/17 domains; 9 sampler, 214 count, 107 numerical dispersion, 14 paired dispersion failure and 959 result-ownership cases |
| Full package checks | Final candidate and unchanged verified original Status OK, including examples/vignettes, `--no-manual` only |
| Actual consumers at 0a | 30/30 stages and complete standard/nuclear Nextflow; exact ordered payloads against unchanged verified original outputs |
| Selected book cells | 14/15; same inherited empty-barcode failure and exit 1 retained |
| Resources | Complete 160 processes at 2ab3f5a plus 48 targeted final 0a processes and reviewed Air source carry; all meaningful timing gates pass; Papalexi cold serial RSS crossing accepted separately |
| Separate costs | Prior eight-million-entry numeric-triplet exception remains bounded; aggregation and dense-copy studies retain their measured pins |
| Maintainability | Four fixed source traces; 7/7 facade-free implementations, 6/7 ordinary-input production R entries; mixed complexity and line-count costs reported |
| Independent source review | Whole 32b source/guide plus separately reviewed D-access and Air deltas; no unresolved P0/P1/material P2 source finding |
| BiocCheck at 0a | Original 0 errors/1 warning/5 notes; candidate 0/1/6; no introduced errors/warnings, one added hygiene note retained |
| Hosted CI at 0a | All six jobs pass in [run 38954047191](https://github.com/Katsevich-Lab/sceptre/actions/runs/38954047191), including Air |

The public API/defaults and 56-slot saved-object schema remain unchanged.
Source/analysis/fit coordinates, storage, QC policy, fitting, score preparation,
cache invalidation/commit, execution and result finalization have explicit
owners. Permutation response traversal/full-run bank and CRT guide traversal/
whole-guide bank remain distinct. Actual consumers use the extracted boundaries;
obsolete routes are retired. [ARCHITECTURE.md](../../ARCHITECTURE.md) explains
concrete future changes without a new extension framework.

The [reviewed behavior ledger](cleanup/behavior-dispositions.md) records named
safety, ownership, import, fit-dimension and plotting corrections.
[User decisions](cleanup/decisions.md) include rejecting repeated file
coordinates, Air formatting and the two separate bounded memory acceptances.
The final Papalexi cold serial peak increases 505 to 667 MiB; the diagnostic
does not establish its cause. The gate remains crossed and further memory
investigation is follow-up work, without restoring input mutation.

Nextflow coverage is the disclosed local shared-arrival configuration, not
uncontrolled ordering/containers/arrays/hosted scheduling. Book coverage is
selected cells. Native probes certify declared ownership boundaries, not real
OOM or every Rcpp allocation path. BiocCheck is not submission approval.
Exact limits and nonblocking follow-ups are in the final report.

Backed CRT/aggregation, SPA, Wilcoxon, accelerated fits, replacement pairwise QC,
Bioconductor data classes, new outcomes, PR233's guide cache and new Python
execution formats remain future work. No new public functionality was required
to finish this internal refactor. Merge/release require a separate instruction.

The detailed execution chronology remains in this file's Git history through
0a9f9d8 and the unchanged unit evidence. Earlier passing counts are historical
evidence, not a substitute for the final report's explicit source/coverage map.
