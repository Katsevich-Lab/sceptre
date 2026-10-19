# Schema 3 architecture

This documents the implemented R object model, not a proposed hierarchy. The
public `sceptre_object` coordinates five persisted components and records
`schema_version = 3L`. Numerical fitting, sampling, scoring, integration and
adjustment remain ordinary functions. `CellView` is a small transient S4 value.
There is no saved Nextflow execution mode or pipeline compatibility claim.

The exact current slot/type inventory is [OBJECT_FIELDS.md](OBJECT_FIELDS.md).
[FIELD_DECISIONS.md](FIELD_DECISIONS.md) accounts for every field in the physical
schema 1 and schema 2 fixtures, including removed fields and migration limits.

## Authority and representations

| Owner | Authoritative payload | Important boundary |
|---|---|---|
| `sceptre_object` | `data`, `spec`, `preparation`, `fits`, `results`, schema version | No duplicate workflow flags or call journal. |
| `ScreenData` | Normalized response/guide counts, source cell covariates, canonical guide metadata, display response names, MOI, imported assignment summaries, content identity | Count-store row names are canonical response IDs; `response_names` is display annotation. Modalities share aligned cells, not features. |
| `AnalysisSpec` | Resolved scientific choices, canonical requested discovery/positive-control pairs, assignment and QC settings | No second `requested` snapshot, budgets, sampled calibration family or worker settings. |
| `PreparationRecord` | Compiled design/pairs/groups, durable guide assignments, selected cells and QC payloads, stage identities | Distinct source and analysis coordinates; assignments exist before QC. |
| `ResponseFitCache` | Compact response coefficients/dispersion plus actual fitting identity | No CRT assignment-model draws or cell-length fitted vectors. |
| `ResultRecord` | Completed entries owning their table, provenance, scope, identity and history; legacy import metadata | A saved table can remain readable while historical. |
| `CellView` | Ordered source/analysis/fit mappings and treatment positions | Built on demand, without copied full count matrices or an eager view per group. |

Canonical requested pair tables in `spec` retain target/response request order
and multiplicity. `preparation@discovery_pairs` and `@positive_control_pairs` are
compiled guide-group expansions. Their `*_with_info` tables additionally retain
nonzero counts and `pass_qc`, including requests that fail QC. Canonical guide
metadata exists once in `data`; preparation stores only its row-aligned
`guide_groups` vector. `guide_group_metadata()` constructs the joined view.

Counts are stored directly as `CountStore = dgRMatrix | BackedCountStore | NULL`;
`NULL` accommodates initialization and recognizable missing legacy stores, not a
usable imported screen. `BackedCountStore` has a typed `CountHandle = externalptr
| NULL` and one descriptor:

```text
BackedCountStore
  descriptor: list(
    path: character[1], content_md5: character[1],
    dimensions: integer[2], dimnames: list(response/guide IDs, cell labels or NULL),
    integer_id: integer[1])
  handle: process-local external pointer or NULL
```

Descriptor dimensions/labels remain available while closed. Numerical access
unwraps a live store and requires reopening. `ondisc` is optional (`Suggests`);
in-memory work neither installs nor loads it. Reopening validates anchored file
content, dimensions, labels and import ID. Verified relocation changes path and
handle without changing scientific identity. The import ID alone proves no
content equivalence. Saving verifies the anchored descriptor and removes live
handles; it does not replace the old content anchor with a new file checksum.
Backing files are immutable during analysis: getters and per-gene loops do not
scan for unsupported external file edits.

## Nested preparation and choice payloads

All indices below are one-based R positions. Ordered identities preserve order;
membership helpers explicitly perform any deduplication or ordering required by
the existing procedure.

```text
ScreenData.import_grna_assignment_info
  max_grna: source-cell-aligned guide labels
  max_grna_frac_umis: source-cell-aligned numeric fractions

AnalysisSpec.grna_assignment_hyperparameters
  maximum: umi_fraction_threshold, min_grna_n_umis_threshold
  thresholding: threshold
  mixture: n_em_rep, pi_guess_range, g_pert_guess_range,
           n_nonzero_cells_cutoff, backup_threshold,
           probability_threshold, formula_object
AnalysisSpec.cellwise_qc_thresholds
  response_n_umis_range, response_n_nonzero_range, p_mito_threshold
AnalysisSpec.additional_cells_to_remove: integer source-cell positions

PreparationRecord.guide_assignments
  <guide ID>: integer source-cell positions
PreparationRecord.source_groups
  grna_group_idxs: named lists of source-cell memberships
  indiv_nt_grna_idxs: named lists of source-cell memberships
PreparationRecord.analysis_groups
  grna_group_idxs: named lists of positions in analysis_cells
  indiv_nt_grna_idxs: named NT memberships
  all_nt_idxs: ordered positions in analysis_cells (NT-control mode only)
PreparationRecord.stages
  design | assignment | cell_qc | pair_qc: list(identity = SHA-256 value)
```

With complement controls, `analysis_groups$indiv_nt_grna_idxs` indexes
`analysis_cells`. With NT controls it indexes `all_nt_idxs`, which itself indexes
`analysis_cells`. `CellView` resolves this distinction into `source_cells`,
`fit_to_analysis`, `trt_fit_pos`, and `n_trt`; counts and design rows use exactly
that order. Complement response fits use all QC cells. NT calibration uses the
ordered NT list; NT discovery/power fits treatment plus NT cells per group and
does not relabel those fits as a global reusable response cache.

## Preparation lifecycle and fitting

The prerequisite chain is import → design → assignment → cell QC → pair QC.
Each preparation stage is committed only after successful computation and
validation. The existence of a matching stage identity, together with valid
prerequisites, establishes completion; there is no saved `complete` Boolean.
`workflow_status(object)` derives public workflow status from these stage
identities and current result records.

Beginning design removes all preparation stage records; beginning assignment
removes assignment/cell-QC/pair-QC records; beginning cell QC removes both QC
records; beginning pair QC removes its record. These transitions mark existing
result entries with `invalidated_by`, retaining their tables and history.
Formula environments are excluded from stable identity values; compiled design
values, dimensions, contrasts/attributes and ordering still participate in the
relevant fingerprints. Failed calls return no partially committed object.

The response cache has a single explicit fitting scope:

```text
ResponseFitCache.entries
  <canonical response ID>: list(fitted_coefs = named numeric vector,
                                theta = positive numeric scalar)
ResponseFitCache.identity
  data: imported content identity
  cells: ordered integer source cells actually fitted
  columns: ordered design column names
  design: SHA-256 of actual selected design values and attributes
  method: "poisson-glm.fit-then-dispersion-v1"
```

Cache lookup is admitted only for a matching fitting identity. Parameter setters,
low-MOI reassignment and retained-cell changes additionally clear incompatible
fits. High-MOI reassignment can retain response fits when their actual fitting
inputs remain unchanged, while tests/results are invalidated. Pair eligibility
can change without changing the cells or response fits. Content replacement
clears preparation and fits and invalidates results.

The algorithm remains Poisson `glm.fit` followed by dispersion estimation at
Poisson means; the method identifier does not claim a negative-binomial refit.
Workers return result/fit proposals. The coordinator validates completed outputs
and commits the resulting value only after success.

## Transient run plans and committed results

A private run plan is a plain list. Calibration preparation returns
`list(object, plan)`; its execution phase takes the plan explicitly. Public
workflow signatures retain their familiar object-in/object-out behavior.

```text
run plan
  analysis: run_calibration_check | run_power_check | run_discovery_analysis
  pairs: complete prepared family, including failed-QC rows
  B1, B2, B3: validated integer budgets
  n_calibration_pairs, calibration_group_size: integer or NULL outside calibration
  identity:
    data, cells, design, assignments,
    choices: side/approximation/control/resampling/integration/reporting settings
    execution: SHA-256 of the plan payload excluding identity
```

Normal plans derive the existing budget policy: B1 = 499; skew-normal B2 = 4999
and B3 = 24999 for permutations or 0 for CRT. Empirical plans use B2 = 0 and the
existing `compute_empirical_budget()` rule. Calibration uses its requested
family size; discovery/power use the maximum of their passing-family counts.
`pair_count(..., passing = FALSE)` counts canonical requested target-response
rows. Passing Bonferroni counts target-response families with at least one
eligible guide; other strategies use their established passing compiled-row
count. A private explicit-budget seam supports numerical fixtures/extensions;
public workflows do not expose arbitrary resampling-budget overrides.

Plans bind actual data/design/ordered cells/assignments, canonical feature IDs,
eligible guide memberships and execution scope. Payload/input changes invalidate
the plan. Calibration choices never mutate working Spec budgets, so a later
calibration cannot overwrite a previous discovery run's realized budget.

Permutation execution shares the full B1+B2+B3 draw bank across the run; CRT
shares a full bank within each guide group. Traversal, native RNGs and draw
consumption remain in the existing numerical drivers.

```text
ResultRecord.records[<public analysis name>]
  table: finalized data frame
  scope: "complete_family"
  identity: frozen run-plan identity
  provenance:
    schema_version, package_version
    spec: value snapshot of resolved choices (formula language, no environment)
    plan: full realized plan and family
    output_amount
    response_fit_method: "poisson-glm.fit-then-dispersion-v1"
  history: list of previous entries, each without another nested history list
  invalidated_by: optional reason
ResultRecord.import_metadata
  migrated_from_schema, provenance, discarded_unverified_fit_count,
  canonical_requests, optional legacy_backend_id
```

Committing requires current QC, a matching plan, full finalized family key
coverage (including legitimate duplicate requests), numeric p-values in [0,1]
or missing, and complete-family adjustment for calibration/discovery. Power
retains its established grouped-result semantics. Noncalibration subset plans
cannot be committed as a complete prepared family. Public analyses always
finalize whole families; there is no live raw-shard result mode.

Currentness requires valid QC, a committed table, `complete_family` scope,
matching current inputs/family, and no `invalidated_by`. Historical validity is
separate: `get_result(..., historical = TRUE)` returns the retained table with
its record attached and warns when it is not current. Plotting, getters and
export use the same current-result guard. Preparing or failing a rerun leaves
the previous committed result intact; successful reruns archive its entry.

Migration inspects physical serialized attributes, not the new registered slot
list. It preserves recognized data/QC payloads, result tables and saved result
history, reconstructs canonical requests only when the saved traversal is
recognizable, and marks old result/fit provenance unknown. It does not certify
old flags, fit reuse, old statistical methods or Nextflow compatibility. See the
explicit retention/removal inventory in FIELD_DECISIONS.md.
