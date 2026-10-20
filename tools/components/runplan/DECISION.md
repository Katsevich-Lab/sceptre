# RunPlan engineering assessment and saved-format decision

A validated list was a reasonable implementation: the old constructor, explicit
input certificate and consume/commit guards already enforced its important
contract. Numerical kernels have no dispatch requirement. S4 is not inherently
more correct or immutable, and the representation change has no claimed speed
benefit. Retaining the list would have been a defensible smaller change.

The requested small `RunPlan` class is also a reasonable domain boundary here:
its eight stable fields are used across preparation, context construction,
execution, result commit and saved provenance. Typed slots and one structural
validity function make that contract inspectable beside the other domain
objects. This is a modest maintainability benefit, not a reason to turn every
nested list into another class. Screen/design/cell/assignment checks still live
at contextual boundaries, and kernels still receive ordinary numerical inputs.

There is one authoritative typed payload. `run_plan_payload()` provides a
canonical temporary value for hashing and legacy interoperation; it is not a
second persisted, mutable settings store. New plans snapshot data.tables that
could otherwise be changed by reference. Ordinary data.frames retain R's
copy-on-modify behavior. Arbitrary slot or reference mutation is not prevented
by S4; consume boundaries verify structural validity and the original payload
certificate before using a plan.

The outer saved schema remains 3. The coordinator/component meanings and the
canonical execution-hash payload remain unchanged; explicit adapters recognize
old schema 3 list plans and response-cache layouts. Known old plans, including
history, convert without recomputing budgets, pairs or input identities. A
real pre-change fixture retains currentness and its original certificate.
Malformed or tampered plans are kept as historical evidence and cannot be
re-signed during conversion. Schema 1/2 unknown provenance stays historical.
