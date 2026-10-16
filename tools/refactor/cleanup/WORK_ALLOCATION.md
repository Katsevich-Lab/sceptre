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
