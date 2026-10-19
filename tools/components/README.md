# Component object model

The public `sceptre_object` coordinates five persistent S4 components and a
schema version. Numerical kernels continue to use ordinary matrices, vectors
and lists. `ScreenData` owns imported counts, aligned canonical cells,
covariates, guide metadata and backend identity. `AnalysisSpec` owns requested
choices; result records freeze the effective settings used by each run.
`PreparationRecord` owns assignment, QC, compiled design and pair eligibility,
with separate preparation stage status. Imported guide metadata remains canonical
in `ScreenData`; derived guide groups belong to `PreparationRecord`. `FitRecord` owns response fits and
identities. `ResultRecord` owns tables and their run provenance. `CellView` is
a temporary, small index mapping used during computation.

Assignments use original source-cell coordinates before QC. Analysis and fit
views validate their ordered mappings rather than copying full count matrices
per pair. Response-fit identity includes data identity, actual design values
and ordered fitting cells; guide probabilities and CRT draws remain transient.
A working-spec change can make a result historical without erasing its table.
Plots, getters and output writers must check current result validity.

Computation returns a proposal before the coordinator commits completion.
Preparing an analysis does not mark it completed; failed computation cannot
commit new results or a completed flag. Worker blocks retain their original
order and sampler behavior. Both association and mixture assignment collection
reject failed workers or lost replies before joining their output.

Schema 2 removes the legacy flat slots rather than maintaining a second copy.
`migrate_sceptre_object()` reads actual serialized legacy attributes and moves
known fields to their component owners. Legacy tables remain readable; unknown
fit and run provenance must not be treated as proof of reusable computations.
Save/read tests cover the real versioned fixture and current round trips.
Backed storage persists descriptors and validates identity when reopening live
handles. Backend location alone is not data content identity.

Nextflow adaptation and compatibility qualification are deferred. Its direct
slot and private-function accesses are not a constraint on this R schema.

The installed ondisc dependency has two preexisting native boundary bugs.
The SCEPTRE pair-count adapter supplies an inert sentinel beyond the feature
row domain to keep an older native helper's final lookahead inside its request
vector; it trims that unused output before restoring request order. Independent
backed count, empty-request and workflow tests assert warning-free results.
The test-only matrix-conversion comparison still calls ondisc's independent
oracle and explicitly accepts its two known pointer-boundary warning messages;
it also accepts zero warnings from a fixed dependency. Other warnings remain
visible and fail that oracle check. These are dependency qualification limits,
not a change to stored-entry counting or statistical algorithms.


Active backing files must remain unchanged for the lifetime of an analysis.
After an in-memory import, the original input files are independent of the
object. Arbitrary external edits to active backing files are unsupported;
getters, plots, reporting and per-gene loops do not monitor or hash the files.
Import anchors file identity, shape, feature/cell names and import ID in
`ScreenData`. Reopen verifies the descriptor and may change the path after
verification. Save preserves the original descriptor and checks it rather than
silently rebaselining claimed original provenance with a fresh hash. Wrong-file
reopen and relocation of unchanged files are covered by regression tests.

In-memory replacement through the internal setters preserves feature and cell
identity, normalizes row access, recomputes count-derived cell summaries and
maximum-guide information, and invalidates preparation/fits/current results.
Replacing backed content or changing feature/cell identities requires a fresh
`import_data()`; legacy mitochondrial feature labels can be unknown, in which
case response replacement also requires reimport. Imported guide metadata stays
canonical while working grouping columns live in preparation.

Requested settings and effective settings are distinct: the original public
parameter arguments are retained, each result freezes its realized budget and
family, and legacy requests are explicitly unknown. Changing the working B3 for
calibration does not alter an earlier run. Power/discovery reconstruct the same
QC budget policy used before this migration. Raw tests are conservatively rerun
on working-spec changes; this implementation does not introduce a raw-test cache.
The response fitter remains Poisson `glm.fit` followed by dispersion estimation,
with a versioned algorithm key that must change when the fitter changes.
