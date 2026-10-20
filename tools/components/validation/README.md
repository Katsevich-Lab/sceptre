# Component-model qualification, 2026-10-19

Implementation commit: `f3c3dc586d0757ee64ba364ec2b80dcbed8939f6`, based on
`97d95da152dbf685e9f4664db7c302844cabe42e`. The source manifest pins package code
and tests; `qualification.json` records versions, artifact identity and counts.
The independent integration review reported no actionable correctness findings
under the documented immutable backing-store contract.

## Newly executed checks

- `R CMD build .`: passed, including vignette generation.
- `R CMD check --no-manual --no-build-vignettes`: **Status: OK**. The built
  vignette and its R code were checked; manual PDF generation was not run.
- Installed-package tests: **4,668 passed, 0 failed, 0 warnings, 0 skipped**.
- BiocCheck 1.47.26 on the exact source tarball: **0 errors, 1 warning, 5 notes**.
  The R-version recommendation and interactive bioc-help check were disabled.
  The warning is the existing absence of Bioconductor infrastructure dependencies.
  Notes concern console output, condition-message construction, function length,
  line width and indentation. These have not been claimed resolved.

Environment: R 4.6.0, macOS arm64, ondisc 1.2.0. Full logs are saved beside this
file. Tests include real low/high-MOI, NT/complement, CRT/permutation workflows;
ordered source/analysis/fit mappings; output levels and failed-QC rows; fit reuse;
assignment/QC/result invalidation; partial preparation and computation failure;
legacy physical-attribute migration; historical results; in-memory replacement
summary refresh; and backed save/reopen/relocation identity.

## Numerical and RNG comparison

Seven newly executed workflows match the unchanged refactor reference at
**tolerance zero with equal realized budgets**, including tables and their
classes, warnings, fit caches, workflow state, and R RNG state. These cover
`low_nt_perm`, `low_complement_perm`, `high_complement_crt`,
`low_complement_crt`, `high_perm`, `low_bonferroni`, and `low_singleton`.
The archived reference library came from the earlier fresh audit; its source
was checked again against the base commit: all 50 R/C++ source files match.
Current native sources have no changes. The comparison runs are new evidence;
the archived reference installation is not a newly built reference.

The raw legacy comparison is also retained. Restoring the original QC budget
policy after calibration intentionally removes a working-budget leak:

| Case | Legacy final B3 | Corrected final B3 | P-values change in this fixture |
|---|---:|---:|---|
| low_complement_perm | 200 | 150 | No |
| low_complement_crt | 400 | 300 | Yes |
| low_singleton | 200 | 250 | Yes |

Other cases match the legacy run directly. In the matched-budget reference,
only B3 is restored externally before power/discovery; reference numerical
functions are unchanged. R RNG states match in all seven raw comparisons too.
The three RDS snapshots, comparison logs, and exact environment-specific probe
scripts in `commands/` preserve this distinction. Original reference outputs
were not overwritten with candidate-derived expectations.

## Scope and remaining limits

The flat saved layout is replaced by schema 2. Public workflows remain familiar;
legacy tables are retained as historical with unknown original provenance, and
unverified legacy fits are not reusable. Arbitrary external edits to active
backing files are unsupported; routine getters/reporting do not hash files.
Save preserves the original descriptor, and reopen checks wrong-file/content
identity while allowing relocation of unchanged files.

Nextflow adaptation and compatibility qualification are deferred. No hosted CI,
Linux/Windows runtime qualification, large-screen benchmark, push, PR update,
merge or release was performed. These small workflows and one real versioned
legacy fixture do not establish compatibility with every historical saved object.
The ondisc private-API boundary and independent oracle warning handling are
explained in the parent design note; no production warnings are suppressed.
