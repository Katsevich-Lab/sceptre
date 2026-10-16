# Cleanup worker ownership and integration

Recorded after independent plan approval, October 16, 2026. This is an execution
allocation under the approved CL units, not a new approval gate. The root
coordinator owns integration, decisions, acceptance, final evidence and PR 234.
Use two focused implementation lanes and an independent senior reviewer; do
not increase simultaneous writers merely to occupy slots.

## Current sequence

| Lane | Worker/worktree | Current ownership and handoff |
| --- | --- | --- |
| Native/count core | rf0_fixture_audit; refactor-ci-oracle/sceptre; codex/refactor-canonical-counts from 489ca7c | CL2: R/storage.R, canonical count/loader native files, relocated group-count block, mitochondrial coordinate ownership, generated bindings and storage/index/ownership tests. CL1 calibration/sampler/inference semantics remain unchanged. Owns test-refactor-indexing.R until explicit handoff. |
| Context/execution | rf2_storage; refactor-rf2-storage/sceptre; codex/refactor-context-execution from 489ca7c | CL3: analysis/context, both drivers/four workers, fitting/score, execution/mixture, three S4 caller migrations, NT-encoder caller and fit-cache source-order contract. Owns obsolete precomputation-file deletion once actual tests migrate; new context tests and the explicit intercept-only assertion repair. Does not touch storage/native/bindings/imports/plots. |
| Independent review | rf0_senior_review, read-only against pinned commits | Review actual source clarity, callers, deletions, indexing/ownership, guards, retained complexity and evidence. Green tests/metrics alone do not establish acceptance. Reviewer owns review artifacts, not implementation. |
| Integration/acceptance | root; integration checkout codex/refactor | Read/review diffs and real contracts, merge coherent tested commits, resolve routine conflicts, own current CI/evidence/decision records and later global tool/name migrations, publish/check CI and qualify final immutable source. |

The count lane released `test-refactor-indexing.R` at
ab4777591cb31906f1c154c1825f0a6bed8c7c6a after its first immutable focused
qualification. The context lane adopts that commit before migrating the real
index/fitting aliases and tests. Count resource acceptance remains open.

Root additionally owns the bounded B1 import repair in `R/import_functs.R` and a
new `test-refactor-import-coordinates.R`: stable complete CSC construction,
actual modality mapping and the explicitly chosen file-coordinate rejection.
Neither implementation worker edits those files. This small independent repair
does not expand the two worker lanes; it requires source review and focused
regressions before acceptance. B6 aggregation and broader import contracts remain
separate gated work.

Root owns the bounded B5 dispersion unit: `src/compute_nb_size_param_functions.cpp`
and the new test-only exception/policy probe and CI step. The executed taxonomy
and independent review support only domain/overflow catches at the existing
whole-estimator try boundary. Controlled MLE/MM-state and unrelated-failure
proof precede acceptance; no production callback, arithmetic or fallback change
is authorized by this allocation. The native retirement lane excludes this file.

Root also owns B5/T7 connection lifetimes in `get_mtx_metadata` and `write_vector`
(`R/import_functs.R`, `R/ondisc_functs.R`) and a dedicated file-lifetime regression
file. Only unconditional helper-owned cleanup is changed; parser formats, normal
file contents and caller-owned connections retain their existing behavior.

Root owns the result adapter/private-alias retirement in `R/result_processing.R`,
the production-only test-utility deletion and its single actual test-helper
caller. Actual pipeline `process_*`/`apply_grouping_to_result` adapters remain;
the unconsumed historical alias and production dimname setters do not. No worker
edits these files concurrently. Root's contract-only import/QC/preparation/source
doc unit is pinned at cacf390; native retirement may adopt those comments before
its subsequent actual-caller name changes.

The conditional B6/T11 experiment uses the existing clean attached worktree
`refactor-execution-cleanup/sceptre`, branch `codex/refactor-vector-aggregation`
from 7df5ad1. It owns only `collapse_grna_matrix`, independent vector-aggregation
regressions and focused resource tooling. After the root's 228 exact
source-evaluated characterization cases, ownership passed explicitly to
rf2_storage for implementation/tests and resource qualification; root does not
edit that worktree concurrently. The branch is not integrated until
exact characterization and the representative allocation/RSS/time gate pass.
The output remains a dense matrix and the backed capability is unchanged.

Root owns the bounded B2/O1 pair-construction unit in
`R/pair_constructor_functs.R`, `src/pair_construction_functs.cpp`, its new
independent regressions and the generated binding for this one native function.
Its distance predicate handles unsorted positions with wide subtraction and
preserves supplied row order; trans exclusions use an owned local copy. No native
sampler/inference names change in this unit. The count lane's current immutable
resource sources are untouched; later native retirement adopts this binding
checkpoint before regenerating the rest.

CL1 at 5ac0fca is integrated as 489ca7c after independent source closure and exact
immutable local qualification (3,631 full checks, 344 focused, 112 native, 12
controlled ownership cases). That prerequisite is complete. The two current
lanes now work on disjoint production files. The context lane prepares its
indexing-test migration as a separate patch until the count lane's coherent
commit releases that file, then rebases/adopts the count change and completes
its own caller/test migration atomically. No broken aliases are accepted.

The count lane coordinates a quiet resource window before large triplet trials;
the context lane may keep editing but suspends competing heavy tests/builds then.
No sampler/inference private renames happen in either current unit: later
retirement jointly updates native definitions, actual R callers, current tests,
ownership probes and tool selectors without aliases. Each implementer documents
functions changed in their owned files. The earlier B1/B6 characterization is
complete and saved with its exact original-source evidence; import/aggregation
implementation follows its own scope and resource gates after this handoff.

## Working rules and handoff

- Workers use isolated active worktrees and named branches from recorded
  integration commits. Preserve unrelated checkpoints and uncommitted work.
- Every task names owned files and exclusions. An overlap/dependency is resolved
  before another writer starts; do not rely on a final merge to discover it.
- Implementers document functions they change against actual callers. No
  separate documentation agent edits the same production files concurrently.
- Each handoff includes full commit, source/build/fixture provenance, focused
  results, known limitations and actual remaining decisions. A successful
  command or queued job is not a passed test report.
- Bindings/actual callers/tests migrate coherently; no new temporary compatibility
  aliases or nullable reconstruction routes are added to split a coupled unit.
- Independent source review plus meaningful contract tests precede acceptance.
  Final broad package/native/consumer/CI/resource gates run on the exact final
  immutable source; historical evidence stays pinned to its original build.
- The coordinator resolves routine engineering choices. Only a genuine unresolved
  scientific, public-interface, duplicate-policy or unexpected blocking conflict
  goes to Eugene; independent authorized work continues meanwhile.
