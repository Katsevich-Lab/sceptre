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
