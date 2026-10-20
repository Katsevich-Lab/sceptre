# Current BiocCheck instance audit and fix scope

Read-only comparison prepared 2026-10-19. Audited current package source: `169c16e4ba00c687c1f03fe37e9d79a11d01cd8f`. HEAD `a3433c71fd339932e4c2d366e5413e13ad8299c9` changes evidence/developer material; its R/man package source is the same. Current source is saved under `source-169c16e/` for stable review. No implementation, Git state, or historical audit artifact was edited.

The coordinator verified the actual submitted source `21f9ea098b69c444a884a49bc254de711968e3b5` against the official source artifact. This report uses `/private/tmp/sceptre-bioconductor-delta-review-20261019/baseline/submitted-source`, the official full `bioc-checks-artifact/00BiocCheck.log`, and the current full `current-00BiocCheck.log`. It does not use PR234 or a pre-refactor approximation as the submission baseline.

## Checker/runtime qualification

The official submission used BiocCheck 1.49.30, Bioconductor 3.24, R 4.6.1/Linuxx86_64 with the official action's `no-check-R-ver=TRUE` and `no-check-bioc-help=TRUE` flags. Its report is0 errors / 1 warning / 3 notes: no Bioconductor dependencies;9 cat/print calls;49 long functions;145 long lines. The current qualified source report uses BiocCheck 1.47.26 on the local Mac and reports0 errors / 1 warning / 7 notes. The coordinator's all-default same-checker baseline additionally records environment/version/subscription checks, including a support-site HTTP timeout; those are distinct from source regressions. Do not add disabling flags to hide current findings.

The focused static predicates here reproduce current69 long functions and473 long lines, and baseline145 long lines. BiocCheck 1.47.26 computes48 baseline long functions, one fewer than the official1.49.30 report. The48-versus49 discrepancy remains explicitly unresolved; it is not grounds to force cleanup of the inherited S4show method without identifying that official49th instance. Upstream current BiocCheck 1.49.33 scanner source was acquired only to inspect implementation; it was not installed or substituted as the submission checker.

## Clear added categories and exact call sites

- `1:3`: `R/distributed_execution.R:290` (`validate_execution_manifest`) and`:343` (`prepare_execution`). Use `seq_len(3L)` or `c(1L,2L,3L)` to retain exactly the accepted output levels.
- Paste nested in a condition signal: `R/execution.R:42`, `announce_parallel_analysis`. Compute the colored pattern before `message()` and pass that value, preserving exact message content.
- Superassignment: `R/perturbation_fits.R:51` (`fit_perturbation_model` warning handler) and `R/schema_migration.R:389,393` (`migrate_preparation_stages/add_stage`). Use a clearly local diagnostic accumulator for warnings while retaining propagation/order, and a helper returning the updated object with explicit rebinding for migration. Do not muffle warnings, remove diagnostics, or bypass validation.
- Missing runnable examples: `man/prepare_execution.Rd`, `man/execute_work_unit.Rd`, `man/collect_execution.Rd`. These are three new exported-object instances. A single shared Rd topic with all three aliases and one real, bounded, installed-package example is appropriate: build/prepare a small valid in-memory screen through public import/setup/assignment/QC; prepare once; execute every unit; collect once; retrieve the completed result. Show that preparation/workers return value artifacts and only collection commits the complete family. `saveRDS/readRDS` can demonstrate transport. Use real public defaults/choices; do not add unsupported private budget knobs, copy source functions, use tests-only fixtures, or insert `dontrun`/`donttest` merely to evade the check. Avoid repeating the full expensive workflow in three independent Rd examples by documenting the lifecycle together.

Exact source, line, column, function identity, baseline mapping, and smallest fix scope for the15 flagged coding call sites are in `flagged-call-instances.tsv`.

## Inherited print/category instances

All9 current cat/print calls have official submission counterparts. `R/output_resources.R:18` relocates the detailed print operation from submitted `R/s4_analysis_functs_2.R:723`. Current S4print cat lines96,98,118,249,300,311,327,340 map to submitted91,93,107,234,266,272,285,295. Retain the resource cleanup fix and detailed textual output. Replacing `print()` with the shorter `show()` output would change behavior and is not required by a no-new-instances policy. The existing no-Bioconductor-dependency warning also remains inherited; do not add an unnecessary dependency to silence it.

## Long function instance scope

`function-lineage-assessment.tsv` lists every one of the69 current flagged function spans, exact file/start/end/physical/coding length, and explicit lineage assessment. `current-all-functions.tsv` and `baseline-all-functions.tsv` include unflagged functions for checking crossings.32 current functions have the same name as a submitted function already over 50 lines. The S4print method and nested plotting helper are clearly inherited despite shifted anonymous labels. Ten renamed scientific/coordinator successors have plausible submitted predecessors, listed with exact old identity/span, but the report marks them as requiring owner confirmation rather than treating name changes as proof. Line-number-suffixed `_anonymous_` labels must never be matched by raw label equality.

The clear new functions/validators and known threshold crossings require either a principled extraction into cohesive, genuinely small helpers or a concrete lineage proof. Do not reduce the metric by minifying bodies, hiding nested implementations in closures, merging statements, removing useful validation, or moving new defects into ignored files. Keep scientific/native algorithms, numeric traversal, RNG/draw sharing, error causality, transactional cache/result admission, and test assertions intact.

Practical extraction seams include distributed manifest certification, unit execution packaging, proposal field/row/scope validation, collector cache merging versus complete-family finalization; schema-one/two copying versus request/history recovery versus identity checks; model-fit entry validation; storage backend counting versus order restoration; and workflow-stage structural checks. A required-API manifest is data: a clearly private constant plus a small accessor can be clearer than a 55-line data-constructor function. Keep constant ownership explicit and do not create mutable public state.

The S4show method is a separate version-sensitive case: local1.47.26 measures50→54lines. It is held unresolved against the official49-function baseline rather than confidently added to the fix list.

## Long line instance scope

`long-line-instance-delta.tsv` lists all 473 current flagged lines with complete text, width, path, line number, whitespace-normalized content, and a multiset-respecting exact-content baseline match where available. It does not offset new instances using removed old lines. The count summary is in `long-line-summary.json`. All unmatched or rewritten lines need wrapping or specific legacy lineage proof; text changes alone are not proof of a new conceptual operation, and an unchanged total would not establish compliance.

Wrap newly added comments/roxygen and call/argument layouts without deleting useful contracts. For generated `R/RcppExports.R`, preserve generated/native-call provenance: use a reproducible formatting step associated with regeneration and verify parsed wrapper/call equivalence; do not change native signatures or algorithms. If a long diagnostic literal must be split, preserve its exact emitted text and avoid introducing another paste-inside-condition finding. The previous AST formatting proof can support pure layout changes, but semantic helper extraction still needs focused meaningful tests plus the required package checks.

## Limits

This is an independent static source-instance audit plus inspection of saved authoritative reports, not a new full BiocCheck/R CMD check or numerical qualification. Existing exact numerical/RNG evidence is not invalidated or expanded by this report. The original and current checkers differ; environment/network diagnostics and the49/48function discrepancy remain separate. A final candidate must be rechecked with the agreed checker/options and an instance-level delta; a category-total comparison is insufficient.
