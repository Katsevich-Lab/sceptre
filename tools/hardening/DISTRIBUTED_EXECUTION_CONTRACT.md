# Distributed execution boundary, protocol version 1

Implementation contract for the separately developed pipeline adapter. No legacy
flat slots are restored. The existing complete-family public workflows remain
the numerical reference and fallback.

- `prepare_execution(sceptre_object, analysis, n_work_units = 1L,
  output_amount = 1L, n_calibration_pairs = NULL,
  calibration_group_size = NULL)` requires completed QC and returns an owned,
  serializable manifest. `analysis` uses canonical `run_*` names. The manifest
  includes a frozen global RunPlan, output level, protocol/runtime/method
  certificates and an ordered `units` list. Calibration pseudogroups are chosen
  here once. Each unit has a stable string `id`, whole response keys for
  permutations or whole guide-group keys for CRT, and original global pair row
  IDs. All eligible rows belong to exactly one unit. No count stores, pointers,
  closures, draw banks or duplicated cell arrays belong in the manifest.
- `execute_work_unit(sceptre_object, execution, unit_id)` validates the manifest
  against the full current object and matching execution runtime. It rebuilds
  full context and local resources, then executes just the certified outer-loop
  keys. Permutation banks use the full assignments, full global cell dimensions,
  global treatment-size bounds and full B1+B2+B3, with the existing fixed native
  seed. CRT keeps each guide's full-budget draws shared across all its responses.
  It returns an owned, certified proposal with explicit global row IDs and raw
  guide-level results plus compact fits. It does not integrate guides, adjust
  across hypotheses, save a ResultRecord, or modify the caller.
- `collect_execution(sceptre_object, execution, proposals)` accepts proposals in
  any arrival order. It validates matching protocol/runtime/plan/context, exact
  expected unit IDs, duplicate/missing/extra rows, unchanged raw row metadata,
  fit scopes and conflicting fits. It restores canonical traversal order, adds
  failed-QC rows from the global plan, integrates guides and applies the complete
  family adjustment once. Only successful collection commits result and fits.
  Failure leaves the input object unchanged.

Returned manifests and proposals are value artifacts. Callers may saveRDS/readRDS
and transport them alongside a separately saved/reopened object. Native/backed
handles and draws are prepared in each worker process, never serialized. The
pipeline must pin one exact package source artifact/image for preparation,
workers and collection. Runtime incompatibility is an error, not an implicit
cache conversion. Certificates guard accidental mismatch/corruption, not hostile
artifact forgery. A protocol or statistical implementation change must update the
protocol certificate and qualify direct versus split/collect behavior.

The build certificate fingerprints actual loaded R function formals/bodies,
package-owned S4 definitions/defaults and registered/effective method code, plus
SHA-256 of the loaded native library bytes. Sorting uses radix order, independent
of locale. It does not trust a supplied environment revision or package version
alone, cache over changed bindings, or load ondisc on memory paths. The exact
source archive/image SHA remains a separate deployment provenance record.

Qualification must include direct versus split permutation/CRT, genuine low-MOI
NT fits, calibration, integration families crossing units, reorder/roundtrip,
missing/duplicate/wrong-context proposals, worker failure and unchanged R RNG.
The adapter remains unqualified against a final candidate until those checks and
its own container/Nextflow checks complete. Backed CRT remains unsupported.
