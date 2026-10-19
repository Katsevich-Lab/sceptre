# Optional ondisc boundary contract

SCEPTRE's in-memory path does not load ondisc. Backed paths inspect actual
capabilities before creating/reopening files or executing native count operations.
A successful check establishes the required R/class/native interface, not a new
numerical or platform qualification. Package version alone is never acceptance.
Checks are cached only for the loaded namespace in the current R process.

## Upstream status checked on 2026-10-19

Authenticated `gh api` reads of the official upstream repository found:

- `main` is [589803edea2a66ad15e7f07186bee16c7739ab8c](https://github.com/timothy-barry/ondisc/commit/589803edea2a66ad15e7f07186bee16c7739ab8c), the pair-iteration bounds fix, with package version **1.3.5**.
- The [GitHub releases](https://github.com/timothy-barry/ondisc/releases) expose only `v1.1` and `v1.2.0`; `v1.2.0` resolves to `5e5e4e55f407ff9025e8ec1c5f84c992b341e2fd`.
- The pinned [NAMESPACE](https://github.com/timothy-barry/ondisc/blob/589803edea2a66ad15e7f07186bee16c7739ab8c/NAMESPACE) exports the public matrix/CellRanger constructors, backing-file reader, example writer and S4 methods. It does **not** export the required pair/group/threshold operations or the constructor accepting a shared import ID.
- The public [matrix constructor](https://github.com/timothy-barry/ondisc/blob/589803edea2a66ad15e7f07186bee16c7739ab8c/R/import_data_from_r_matrices.R) forwards to an internal constructor whose default `integer_id` is zero. Its exported formals do not let SCEPTRE supply the shared ID.

Thus the inspected official upstream API still leaves a support/distribution gap:
there is no discovered supported released replacement for the private contract
below. This is scoped to the official GitHub source/releases inspected, not a
claim about a new release, a Bioconductor ruling, or an upstream commitment.
No upstream contacts or library upgrades were performed.

The recommended source remains the exact README/CI pin above. Earlier project
qualification identified a same-version 1.3.5 source without the pair bounds
fix, so `Version: 1.3.5` cannot replace a source pin. The local installation used
for this hardening test is separately **ondisc 1.2.0**, R 4.6.0/macOS arm64, with
no `RemoteSha` in DESCRIPTION. It passed the interface gate and the existing
backed semantic regressions; its version does not prove it is that pinned source.
No new installed-package qualification of upstream main or another release is
claimed here.

## Required operations

`R/ondisc_capabilities.R` is the executable interface manifest. It verifies the
following function argument lists, the exported reader, the required `odm`
slots, and registered native arities. Private calls remain confined to storage.

| Operation | Required entry | Required semantics |
|---|---|---|
| Create counts | private `create_odm_from_r_matrix_internal(mat, file_to_write, chunk_size, compression_level, integer_id)` | Preserve supplied shared integer ID; row-accessible counts, feature order and integer-conversion behavior. The response and guide stores share the caller's import ID. |
| Reopen | exported `initialize_odm_from_backing_file(odm_file)` | Create a process-local handle and expose saved dimensions, feature IDs and import ID. Both file IDs must match before use. |
| Row read | `odm` row `[` method, registered `load_row_cpp` | Full gene row in original column order, followed by SCEPTRE selection in requested source order; repeats remain repeats. |
| Assignment threshold | private `threshold_count_matrix_cpp` | Original one-based stored positions; ondisc coerces count/threshold to integers. Stored zeros can pass threshold zero; unstored zeros never appear. |
| Pair counts | private `compute_n_ok_pairs_ondisc` | Count stored positions in treatment/control groups for requested rows; selected-cell and NT coordinate rules remain unchanged. |
| Pair + NT summary | private `compute_nt_nonzero_matrix_and_n_ok_pairs_ondisc` | Same pair counts plus NT-group-by-gene counts and totals, including supplied group overlaps. |
| Group counts | private `compute_n_trt_cells_matrix_ondisc` | Preserve group order, requested response order/repeats, dimensions and supplied membership multiplicity. |

Required raw `odm` slots are `h5_file: character`, `dimension: integer`,
`feature_ids: character`, `ptr: externalptr`, and `integer_id: integer`.
The typed SCEPTRE `BackedCountStore` remains independent of ondisc at package
load: its descriptor is the sole saved authority and its optional handle is
process local. Public count getters reconstruct the existing raw `odm` view;
numerical kernels and per-gene count paths are unchanged.

SCEPTRE owns copies of group-coordinate vectors before calling native routines
that decrement them in place. It sorts pair requests for reading, restores
original request order afterward, and appends an inert `n_genes + 1` request to
keep legacy native lookahead within bounds. The sentinel is never read as a gene,
never enlarges the statistical family, and its unused output is discarded.
These safeguards remain even when a compatible installation includes the bounds
fix. This change introduces no alternate count kernel, RNG consumption or CRT
support for ondisc.

## Saved-object and content lifecycle

The backing descriptor contains absolute path, file MD5, dimensions, feature IDs
and shared integer ID. The ID establishes pairing, not immutable content.
Saving/reopening verifies the anchored checksum; verified path relocation does
not invalidate compatible fits/results. Public backed saving stores closed
handles. `readRDS()` alone permits descriptor and saved-result inspection but
cannot execute count-dependent work; the public backed reader opens current
process handles. See the read/write help for runnable save/read examples.

Active backing files are immutable. Getters and per-pair loops do not scan full
files to detect unsupported external edits. Changed count content/import
summaries require reimport. A capability failure explains the installed version,
missing or changed operations, qualified source pin, restart requirement and
availability of the independent memory path; it occurs before creating count
files. Checking capabilities cannot prove that arbitrary new native code keeps
these semantic promises: qualification against these tests remains necessary.

## Regression evidence and remaining qualification

- `test-ondisc-capabilities.R`: compatible API independent of version text;
  missing operation, incompatible class, missing exported reader, changed
  reader formals/native arity; actionable absent/incompatible diagnostics;
  failure before files are created; installed API acceptance.
- `test-refactor-storage.R`: stored-zero/fraction differences, row/column
  ordering, owned native coordinates, sentinel behavior, group counts and
  shared-ID save/reopen checks.
- `test-schema-migration.R`: typed closed stores, relocation/tamper rejection,
  metadata-only descriptor methods, legacy migration and saved results.
- `test-refactor-contracts.R`: actual backed import/QC/save/reopen/discovery
  agreement with sparse workflows and wrong-file-ID rejection.

Incompatible-install tests use injected metadata snapshots, not suspect native
code or replacement packages in the ordinary library. The cold-process probe
loads the package, migrates an in-memory fixture, and inspects a closed typed
store without loading ondisc. No broad platform, R-floor, backend CRT or new
upstream-release qualification follows from these local checks.
