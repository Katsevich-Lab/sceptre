# Actual Bioconductor submission baseline and required delta scope

Status: the actual submission source and original official diagnostic report are established. The hardening source 169c16e is **not** within the requested no-new-diagnostics boundary. Final implementation changes and an exact committed candidate are pending independent recheck.

## Authority and exact submitted artifact

The authoritative submission is [Bioconductor/BiocContributions issue #176](https://github.com/Bioconductor/BiocContributions/issues/176), submitted August 31, 2026. Its first-build comment links the original [run #33438717492](https://github.com/r-universe/biocstaging/actions/runs/33438717492) and [Bioconductor-check job #99642866650](https://github.com/r-universe/biocstaging/actions/runs/33438717492/job/99642866650).

- Actual submitted source: `21f9ea098b69c444a884a49bc254de711968e3b5` in BiocStaging/sceptre, branch devel. The pinned run tree at `77f084cc9ed489918bad9d413821dac11a037a37` records this exact submodule SHA, and both original job inputs/source-build CommitRef agree. Authenticated current staging history still points to21f9ea.
- Original official source archive SHA256: `5300b3e3adfbb73a9a2173384ea3ec7e57b8cad8612f3fdc8116d5db5156eaf4`.
- Retained source-artifact ZIP ID 9775542189 SHA256: `59fd0f2758857f648f9eb74140188b6e8fc77b75fa75bdf2279bf54f19a267b8`.
- Retained BiocCheck-report ZIP ID 9775660207 SHA256: `96d7662cd2e75abb4f316f503a8de2c59fc3f02a4d4409ca19e48bb11bf79d05`.
- Both original artifact ZIP digests were independently verified. All 97 compared packaged production/documentation files match the submitted Git source byte for byte (`baseline/submitted-source-equivalence.json`).
- Original checker: **BiocCheck 1.49.30**, Bioconductor 3.24, package0.99.0, Linux x86_64 / R 4.6.1. Original container reference was ghcr.io/r-universe-org/base-image:release, resolved digest `sha256:9e2813ce8cf2db1e107b2b195481f14b5c08d828c5d3692707e280aaef22cd14`.

The later refactor PR234 and its0E/1W/5N qualification are not the submitted baseline. The local submission plan is useful context but its local five-note header is not the original official result. The separately referenced cover-letter Markdown file is absent at the supplied location; the original official report/source are nevertheless available and sufficient to establish this baseline.

## What the original official report actually allowed

The original official report is **0 errors, 1 warning, 3 notes**:

| Diagnostic | Original flagged scope |
| --- | --- |
| No Bioconductor dependencies | One warning; the public September 1 explanation discusses the specialized perturb-seq representation. |
| cat/print outside show | Nine calls: one summary-serialization print in s4_analysis_functs_2.R and eight cat calls in s4_helpers.R. |
| Function length |49 functions over 50 lines in the official1.49.30 report; only the longest 5 names are printed. |
| Long lines |145 lines over 80 characters, with the full flagged-location list in00BiocCheck.log. |

Preserved originals: `baseline/bioc-checks-artifact/00BiocCheck.log`, `output.log`, original Actions logs plus ANSI-stripped reading copies, original source tarball and pinned wrapper scripts. The [public September 1 comment](https://github.com/Bioconductor/BiocContributions/issues/176#issuecomment-5486819338) corroborates the three source-note categories and 145 long lines; it is not a substitute for the original artifact.

## Checker and environment distinctions

The pinned official wrapper at r-universe-org/actions commit `e859fe2adc2553a5c8b0b2ec6044215d60d3d9c7` historically passed `no-check-R-ver=TRUE` and `no-check-bioc-help=TRUE`. That exact original script is preserved; these are historical report settings, not newly adopted suppressions.

This audit also ran the **original official source archive through installed BiocCheck 1.47.26 under R 4.6.0 with all checks enabled and no disabling flags**. Its terminal result is **1 error, 1 warning, 5 notes**. The extra error is a support-site HTTP timeout; the two extra notes concern the unchanged R>=4.1 declaration and mailing-list status requiring administrator credentials. They are environment/checker-setting differences, not hardening-source regressions. This run is not reported as clean. Its source categories remain the same three, with 48 long functions and 145 long lines. The49-versus 48 function count is explicitly version-sensitive and must not be silently normalized.

Files: `baseline/same-checker-baseline.log`, `baseline/same-checker-1.47.26/diagnostics.json`, `result.rds`, `session-info.txt` and `invocation.txt`. Final candidate source must be compared under the same checker/options, retaining any environmental failure rather than disabling checks to hide it.

## Required implementation repair scope

The pre-fix hardening source is 169c16e, with packaged R/man files unchanged at initial owner HEAD a3433c7. Its existing qualified1.47.26 report is0E / 1W / 7N under the historical two exclusions. Four new source-note categories are clear:

1. Two `1:3` output-level membership expressions in `R/distributed_execution.R` (lines 290, 343).
2. `paste` inside a condition signal in `R/execution.R` (line 42).
3. Three `<<-` sites: `R/perturbation_fits.R` line 51 and `R/schema_migration.R` lines 389, 393.
4. Missing runnable examples on prepare_execution.Rd, execute_work_unit.Rd and collect_execution.Rd. Examples must actually run installed prepare→execute-each→collect on a small real dataset; donttest/dontrun/comment-only examples would not meet this requirement.

An inherited category does **not** exempt fresh violations:

- All nine current cat/print operations have direct submitted counterparts. The current output_resources.R print is the relocated summary-serialization print; the eight current s4_helpers calls map to the original display method. They are inherited instances, so replacing detailed output with a shortened show() or logging would be unnecessary and could change behavior.
- Current1.47.26 flags 69 long functions versus 48 for the submitted archive under that same version. Every new or newly over-threshold function needs a meaningful decomposition or documented source lineage to an already-flagged submitted function. Anonymous functions must be matched by enclosing method/helper, never by line-suffixed checker labels. Counts alone cannot establish equivalence.
- Current formatting flags 473 long lines. 110 match whitespace-normalized original long-line contents; 363 are new or rewritten candidates requiring wrapping or specific inherited-lineage evidence. A final total<=145 would be insufficient if it still contained new flagged lines.

The machine-readable source/function/line inventory lives in `current-instance-audit/`; its lineage report distinguishes verified inheritance from unresolved/new instances. Line numbers here identify the frozen pre-fix source and are not candidate locations.

## Final acceptance evidence required

Identify the exact committed candidate and archive hash. Preserve the full diagnostic reports, checker versions, options and environments for submitted and final source. Demonstrate zero newly introduced flagged instances in addition to no new source-note categories; do not use suppressions, meaningless dependencies, dense code packing or unexecuted examples. Independently run the actual installed examples and inspect changes for statistical/native behavior, RNG/draw order, complete-family correction, lifecycle/migration/cache identities and public-workflow preservation. Existing source/runtime/image qualification is historical until it binds the new candidate; the coordinator owns any required pipeline/image refresh.

Retrieval scope: one exact live Compass project note; authenticated issue #176 and its 9 comments; pinned official run/job metadata and full logs; two original retained artifact ZIPs; pinned checker-wrapper code; current staging history; and the supplied local submission tracker/HTML. This is not a claim to have reviewed all project knowledge. No implementation, other repository, original submission evidence, external maintainer communication, push, merge or release was modified/performed.
