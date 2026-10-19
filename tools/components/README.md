# Component object model

The final schema is version 3. `sceptre_object` coordinates five persistent
components: imported `ScreenData`, requested `AnalysisSpec`, compiled
`PreparationRecord`, compact `ResponseFitCache`, and entry-owned `ResultRecord`.
`CellView` is a temporary computation value. A typed per-run `RunPlan` is created
for execution and retained with completed result provenance.

- [Exact per-object fields and declared types](schema3/OBJECT_FIELDS.md), generated
  from the current S4 class definitions with exhaustive description checks.
- [Architecture and nested payload contracts](schema3/ARCHITECTURE.md).
- [Exhaustive legacy-field ownership decisions](schema3/FIELD_DECISIONS.md).
- [Schema 3 qualification](schema3/VALIDATION.md).
- [Typed RunPlan checkpoint qualification](runplan/validation/README.md).
- [Earlier schema 2 qualification](validation/README.md), preserved as historical evidence.

Canonical requested target-response pairs stay in `AnalysisSpec`. Expanded guide
families, design, assignment and QC payloads belong to preparation. Preparation
stage presence plus its validated input identity establishes completion; status
and plot selection are derived. Cheap counters and last-called flags are absent.
Run budgets, realized calibration families and execution choices are transient;
a completed result freezes its effective plan and requested-spec snapshot.
Each result entry owns its table, identity, scope, provenance and history.

Count slots use the `CountStore` union of `dgRMatrix`, `BackedCountStore` and
`NULL`. The small backed adapter owns one descriptor and an optional live
external pointer. It does not name any optional class in its type definition.
Package loading, in-memory workflows and closed descriptor inspection do not
load `ondisc`; only backed numerical access requires it. Save closes the handle
while retaining the anchored descriptor. Reopen verifies content identity and
allows relocation of unchanged files.

Active backing files must remain unchanged during an analysis. Getters, plots,
reporting and per-gene loops do not hash them. In-memory replacement preserves
feature/cell identities, refreshes import summaries and invalidates dependent
preparation, fits and results. Replacing backed content requires reimport.
Schema 1/2 migration reads physical saved attributes only at the boundary;
unknown results remain historical and unverified fits are discarded.

The existing ondisc pair-count boundary uses an inert sentinel to accommodate
an older native helper's final lookahead. Production tests require warning-free
results. The independent matrix-conversion oracle accepts only its two known
pointer-boundary warning messages, or zero warnings from a fixed dependency.
Native statistical kernels and samplers are unchanged. Nextflow adaptation,
large-screen benchmarking and cross-platform qualification remain separate work.
