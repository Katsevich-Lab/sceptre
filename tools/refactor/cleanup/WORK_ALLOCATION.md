# Cleanup worker ownership and integration

Recorded after independent plan approval, October 16, 2026. This is an execution
allocation under the approved CL units, not a new approval gate. The root
coordinator owns integration, decisions, acceptance, final evidence and PR 234.
Use two focused implementation lanes and an independent senior reviewer; do
not increase simultaneous writers merely to occupy slots.

## Current sequence

| Lane | Worker/worktree | Current ownership and handoff |
| --- | --- | --- |
| Native/core | rf0_fixture_audit; refactor-ci-oracle/sceptre; branch codex/refactor-calibration-safety from 6a43fa5 | CL1 only first: calibration validation/derived budgets/permutation preparation, sampler/calibration native safety, actual inference boundary/stage/tail safety, focused new regressions. Count/group kernel restructuring is excluded from this initial unit. No broad private renames yet. |
| Independent contract preparation | rf2_storage; refactor-rf2-storage/sceptre | B1 import/duplicate and B6 in-memory aggregation characterization against immutable original libraries; no production changes while CL1 is pending. Small safe probes only, with no invalid native layouts or competing heavy resource jobs. The repeated-file-coordinate choice is recorded in decisions.md. |
| Independent review | rf0_senior_review, read-only against pinned commits | Review actual source clarity, callers, deletions, indexing/ownership, guards, retained complexity and evidence. Green tests/metrics alone do not establish acceptance. The reviewer owns review artifacts, not implementation. |
| Integration/acceptance | root; integration checkout codex/refactor | Read/review each diff, check real contracts, merge coherent tested commits, resolve routine conflicts, regenerate joint bindings where needed, publish/check CI and qualify the final immutable source. |

CL1 must pass implementation review and focused compatibility/safety gates before
broad CL2–CL6 restructuring. After that handoff, the native/core lane can own
canonical counts/index ownership and native retirement. The second implementation
lane can own the tightly coupled context/driver/worker migration or an independent
import/aggregation unit, depending on finished prerequisites. The coordinator
records concrete file ownership at each dispatch; one lane owns each coupled
migration atomically. No simultaneous edits to shared context/storage/fitting,
Rcpp bindings or facade callers.

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
