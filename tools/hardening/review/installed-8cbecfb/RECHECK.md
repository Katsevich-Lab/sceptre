# Independent focused recheck at 8cbecfb

**Both open R-boundary findings are verified fixed at `8cbecfb88afdfbb3b45a991f63def7bcef8c8848`.** This is a disposition for that exact source, not final cross-repository signoff. No implementation/Git edits were made.

I exported the exact commit with `git archive`, installed it from clean source with `R CMD INSTALL --preclean` into a new isolated library, and ran the original self-contained reproducer in **installed mode**, followed by12 focused checks. This used local macOS arm64 / R4.6.0 and a freshly compiled native library; it did not overlay the old namespace or reuse its compiled package.

| Finding | Disposition and independent evidence |
|---|---|
| Conflicting warm perturbation fits accepted | **Fixed.** A genuine matching CRT proposal collects and remains current. The original re-certified changed coefficient now rejects without changing the caller. Changed diagnostics and an attempted NULL removal of a valid same-input warm fit also reject atomically. A genuinely stale per-entry treatment identity misses and is successfully replaced, so the conflict guard does not prevent legitimate refresh. |
| Full ordered cells repeated in proposals | **Fixed.** The same returned proposal graph walk now finds zero full ordered-cell paths in either shared-fit CRT proposal. The global manifest retains its one authoritative cell vector. Full response/perturbation scope identities are reconstructed correctly in the committed caches. |

Post-collection attribute mutations through caller-retained response and perturbation proposal coefficient vectors remain isolated. All12 extended checks passed. The source diff explicitly implements the NULL policy—valid matching warm entries cannot be removed—and its tests cover coefficient, diagnostic, removal and stale replacement cases.

The three repaired manually constructed fixtures now explicitly supply valid `multiple_testing_alpha=0.1` and `multiple_testing_method='BH'`. The inspected diff does not weaken the initialized-spec validator or change those tests' coordinate/legacy-behavior expectations.

## Exact artifacts

- Source archive SHA256: `2c84adca90818abb8cc55ef049a4059d2a4214e564a6d0d4e0bd2089aa50210a`.
- Distributed R source SHA256: `04780067a6997855c93ea1da576bdd7e0900fa12a93e8c57f521c6c7b7615014`.
- Fresh native library SHA256: `a56a9cb23f44c6a793d35b030d08d7e7817f31d41b6e12056c997e81f8e8423a`.
- [Artifact manifest](/private/tmp/sceptre-hardening-review-20261019/installed-8cbecfb/artifact-manifest.json): source and all installed-file hashes.
- [Machine-readable verification](/private/tmp/sceptre-hardening-review-20261019/installed-8cbecfb/verification.json): exact loaded runtime/build certificate and12 results.
- [Install log](/private/tmp/sceptre-hardening-review-20261019/installed-8cbecfb/install.log), [original repro recheck](/private/tmp/sceptre-hardening-review-20261019/installed-8cbecfb/recheck.log), [extended recheck](/private/tmp/sceptre-hardening-review-20261019/installed-8cbecfb/recheck-extended.log).
- [Original runner](/private/tmp/sceptre-hardening-review-20261019/open-findings-reproducer/reproduce.R) and [extended checks](/private/tmp/sceptre-hardening-review-20261019/open-findings-reproducer/extended-recheck-8cbecfb.R).

The fresh library is `/private/tmp/sceptre-hardening-review-20261019/installed-8cbecfb/lib`. These results are not a broad suite, R CMD check, R4.4 qualification, performance benchmark, container run or actual Nextflow qualification. The earlier27-case numerical split/collect comparison is separate evidence for its earlier source snapshot. Final integration must bind the later frozen package/pipeline revisions, image and runtime, and recheck any intervening relevant source changes.

The independent pipeline finding remains open at pipeline `5e257b5`: the run signature does not bind Nextflow implementation changes. Its [separate recheck report](/private/tmp/sceptre-hardening-review-20261019/pipeline-review/RECHECK_5e257b5.md) records the passing ordinary cases and the remaining actual module-change/resume failure. That issue is not closed by these R-package fixes.
