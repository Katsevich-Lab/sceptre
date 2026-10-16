# Cleanup decisions after independent plan approval

The reviewed target is pinned at 08e0260 and its approval in
senior-review-amendment.md. This register records subsequent explicit user
choices; it does not erase the original findings or qualifications.

## D2026-10-16-01: repeated imported coordinates

Eugene selected: **Reject repeated coordinates with a clear import error**.
The question concerned repeated `(cell, feature)` entries in imported count
files, for which preserving separate entries cannot produce a valid CSC matrix.
The alternative of summing entries was presented and was not selected.

This resolves B1/T1's duplicate-policy stop. The file-import boundary must detect
repeated coordinates and fail informatively before constructing/using compressed
sparse slots. Valid unique coordinates are sorted stably and modality IDs mapped
explicitly, with values, explicit zeros, dimensions and IDs preserved. No silent
coalescing, summation or last-entry selection is introduced.

This decision does not change valid triplet-matrix duplicate handling or repeated
guide-metadata multiplicity in the existing in-memory vector aggregation path.
Those are separate contracts being characterized under B6/T11. Add independent
file-import rejection and unchanged ordinary-input regressions; never run the
invalid old CSC layout through downstream native kernels.

Implementation boundary: B1 repairs SCEPTRE's owned in-memory Cell Ranger and
Parse constructors. The `use_ondisc=TRUE` Cell Ranger route delegates file
parsing directly to ondisc; this unit neither changes nor qualifies that
external parser's duplicate/modality behavior. Independent source review
requires this limit to remain explicit. A full-file in-memory preflight would
violate the backed route's resource boundary and is not introduced here.

## D2026-10-16-02: measured numeric-triplet ownership cost accepted

Eugene explicitly replied **“I accept”** to the completed revised study's
bounded correctness-cost disposition. The disclosed experiment compares
original 21f9ea098b69c444a884a49bc254de711968e3b5 with candidate
909aa3ad76ee21e7f9ab83a63c645bc0fb6a43c7. At eight million stored entries,
numeric row/column conversion adds **121.938/123.141 MiB** whole-child peak RSS
and crosses the generic memory gate; logical adds **61.891/58.188 MiB** and
passes. Exact outputs, caller nonmutation and all meaningful timing gates pass.

**Disposition: accepted bounded correctness-cost exception for those measured
numeric-triplet cases. The generic gates remain crossed, not passed.** Preserve
all initial/revised measurements. Stop optimization or reruns whose sole purpose
is making these unchanged numeric cases green. No additional user confirmation
is needed for the disclosed cost.

The mechanism owns borrowed coordinates/values before in-place sorting, avoiding
caller mutation. Constructor-proven fresh logical values are already reused.
This is a source/fixture/platform-specific acceptance, not a universal bytes-per-
entry allowance, proof of an algorithmic minimum, larger-input memory bound or
waiver of unrelated regressions. Later source changes materially affecting this
cost/behavior require delta assessment. Correctness, nonmutation, index safety,
evidence/source consistency, independent mechanism review, the Linux aggregation
precision repair and final integrated qualification remain required.

Authoritative report SHA256:
`9dfd92d6cf625b8889da8931f5f811488fef32b7a8a86cda98c9564db7db2868`.
[Losslessly retained full study and separate dense-cost evidence](../results/conversion-cost-909aa3a/README.md).
The original proposal and study's pre-acceptance wording are preserved; this
entry supplies the subsequent human disposition without changing their bytes.

## D2026-10-16-03: Air formatting authorized

Eugene explicitly selected Air as the single formatter and authorized applying
it. Retain the existing four-space indentation and preferred width of 80,
pin Air 0.9.0 in the formatting CI check, and wrap overlong prose comments and
documentation separately where useful. Width is a preference, not a guarantee;
BiocCheck remains a distinct checker.

This supersedes the earlier nonblocking style deferral only for this formatting
unit. Scope is handwritten `R/`, formatting configuration/contributor guidance,
a scoped CI check and this authorization record. Preserve literal values,
formals, help meaning and mathematical behavior. Generated bindings, immutable
fixtures, reference/tooling/results and tests are not recursively formatted.
Intentional long code/example/URL/string lines may remain; do not rewrite
executable diagnostics or unrelated code to meet a line limit. Audit parsed
expressions with narrowly normalized single-expression braces and inspect
meaningful deltas. The combined new source pin still requires qualification;
earlier resource measurements remain evidence for their measured pin.
