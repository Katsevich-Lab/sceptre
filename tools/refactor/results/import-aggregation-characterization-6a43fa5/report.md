# B1 and B6/T11 bounded characterization

Status: **read-only characterization completed; no production changes, resource
qualification, new backed aggregation or general bug-freedom claim**.

Plan authority: `6a43fa5b2d32109a67e42138d3d845a6e5a2805a`, including the approved
behavior dispositions and senior amendment. Executed original source:
`21f9ea098b69c444a884a49bc254de711968e3b5`; installed library:
`/private/tmp/sceptre-resource-reference-release-library`. The historical
28df390 candidate library was not used or relabeled as final qualification.

Run `sh /private/tmp/sceptre-import-aggregation-characterization-6a43fa5/run.sh` to reproduce the small probes.
`characterize.R` and `rounding_sweep.R` use separate independently copied inputs
for original and scratch calculations. The report pins source files, installed
function formals/bodies, dependencies, fixture hashes and DLL hashes. Exact
numeric values are recorded as decimal and hexadecimal doubles in report.json.

## Provenance and scope

- R version 4.6.0 (2026-04-24), `aarch64-apple-darwin23`; source ondisc 1.3.5.
- Five installed original function formals/bodies matched their pinned Git source.
- Original sceptre DLL SHA256: `3caa868fd34acf0bb04b6ab3ccc418718e913c8ab76aa4bf53facfbbdda8c37e`; it matches the
  verified clean release-build manifest, whose complete compiler/source evidence
  is embedded and hash-pinned in report.json.
- Matrix, data.table, ondisc, dplyr, digest and jsonlite versions/paths are recorded.
  Matrix/data.table/ondisc DLL hashes are also recorded. Thread environment is
  fixed to one; data.table reports one thread. Public Cell Ranger source itself
  requests nThread=2 in its tiny fread call; this is characterization, not a timing
  comparison or a claim to have changed that reader.
- All probes are small; no heavy concurrent profiling, large allocation, package
  rebuild or unsafe native baseline example was run. The artifact directory is
  under 1 MiB at this checkpoint.

## B1: exact conflict and recorded decision

A raw `dgCMatrix` with `i=c(0,0)`, `p=c(0,2,2)`, `x=c(2,3)`, `Dim=c(2,2)`
retains both duplicate entries but safe `Matrix::validObject(test=TRUE)` returns
**"'i' slot is not increasing within columns"**. The corresponding sorted unique
`i=c(0,1)` is valid. Unsorted unique rows within a column also fail that check.
Stable sorting alone cannot preserve duplicated stored entries and make this
CSC representation valid. These malformed layouts were checked only for validity
and their raw slots; they were never converted, combined or passed to downstream
statistical/native kernels.

**The policy choice is resolved.** Parent relayed Eugene's explicit decision:
**reject repeated file (cell, feature) coordinates with a clear import error**;
do not coalesce or sum them. That authority is specific to file-coordinate B1.
It does not change valid dgT duplicate aggregation or repeated metadata-row
multiplicity. No new error is implemented by this artifact. A repeated coordinate
within one file differs from the same feature in a different cell or independent
batch, which remains valid.

Validity is necessary but not sufficient for import completeness. For input cell
indices `c(1,0,1)` and unique feature indices `c(0,1,2)`, source reasoning gives the
old pointer `c(0,0,1)`. Matrix accepts raw slot capacity greater than the last
pointer, so validObject returns TRUE although only one of three supplied entries
is visible. No downstream operation was run on that layout. A repair must assert
all stored entries/pointers against independent expectations, not merely validity.

The owned, stable intended ordering example has input rows
`(cell,feature,x) = (1,1,11),(0,1,1),(1,0,10),(0,0,0)`. The independent order is
`4,2,3,1`, producing `p=c(0,2,4)`, `i=c(0,1,0,1)` and `x=c(0,1,10,11)`.
The explicit zero remains stored. This is an intended scratch calculation, not
an executed candidate import implementation.

Interleaved feature IDs show the separate modality bug without executing it:
genes at global rows `1,3,5` become old rows `0,2,4`, outside a three-row matrix;
guides at `2,4` become `0,2`, outside a two-row matrix. Actual matching to the
modality's ordered global-ID vector gives gene `0,1,2` and guide `0,1`. The old
out-of-range layouts were subjected only to safe validity checks. Preserve
feature-file row order, not alphabetical feature-ID order.

A safe ordinary two-directory public Cell Ranger import succeeds. Its returned
response and guide matrices are `dgRMatrix`; rows follow feature-file order;
columns/covariate rows are exactly `b1_cell-z,b1_cell-a,b2_cell-z,b2_cell-a`, with
batch `b1,b1,b2,b2` following supplied-directory order. The final matrices retain
the ordinary input's stored zeros/fractions. Grouped unique create_csc_matrix
also has its actual returned class/pointers/indices/values recorded. No repeated
file-coordinate or interleaved-modality input was sent through the old unsafe path.

## B6/T11: observed current contract

There are **28/28 exact original-versus-scratch output matches**, with every
original and scratch caller input unchanged. This covers numeric/integer/logical
dense matrices; valid numeric and logical C/R/T sparse matrices; unsorted valid
triplets; stored zero/FALSE entries; valid numeric/logical duplicate triplets;
repeated and cross-vector metadata; one cell/vector; helper-only zero cells and
zero vectors; extra metadata columns; factors and missing vector labels; and
fractional, large-rounding, cancellation/overflow helper cases. Helper-only edge
cases do not assert a new supported public workflow or add validation.

For the hand-specified normal fixture, source guide order is
`g3,g1,g4,g2,unused`, metadata first-seen vector order is `vZ,vA,vB`, and cell
order is `cell-z,cell-a,cell-m`. Original helper output is a dense numeric
`matrix`/`array`:

| Vector | cell-z | cell-a | cell-m |
| --- | ---: | ---: | ---: |
| vZ | 4 | 0.1 | 4.5 |
| vA | 6.5 | 4 | 2.5 |
| vB | 0 | 3 | 0.125 |

`vA` includes `g1` twice because its metadata row is repeated, and `g2` contributes
independently to both vZ and vA. The unused guide contributes nowhere. Do not
`distinct()` guide members or collapse metadata before computing these sums.

The separate safe duplicate example is exact: valid `dgTMatrix` stores g1's
values 2 and 3 at the same cell; the existing Matrix densification makes g1=5.
Two g1 metadata rows then yield vector v1=10; g2 contributes v2=7. Output order is
v2,v1 with cells c2,c1. Valid logical triplet duplicates use the existing logical
value semantics, then metadata repeats contribute independently; do not count
stored TRUE duplicates as new numeric counts.

`collapse_grna_target_data_frame` drops guide ID, renames vector ID to guide ID
and takes full-row distinctness in first-seen order. Exact duplicate metadata
rows disappear in this **returned metadata**, even though they counted twice in
aggregation. Different extra-column annotations can leave multiple metadata rows
with the same vector ID. The matrix still has one row per first-seen vector.
These are distinct multiplicities and must not be conflated.

Original helper output remains dense, even for sparse input. Public in-memory
import subsequently converts it to `dgRMatrix`; four public numeric/logical sparse
vector imports succeeded and record exact final matrix and metadata state.
Aggregation's dense zeros do not retain an input sparse-position identity;
subsequent dense-to-sparse conversion drops zeros. This differs from the
nonaggregating file-import stored-zero path and is existing behavior.

The helper retains its simplification/type details: zero-vector output is 0 by
n_cells with no dimnames, a one-vector result remains a matrix, factor vector IDs
can yield no row names, and factor guide indices retain the original indexing
behavior. Missing vector IDs are filtered according to the original dplyr
selector. A casual replacement with base logical indexing would treat NA masks
differently. The scratch proof therefore keeps the original dplyr filter/pull
selector and t(vapply()) output construction unchanged; it does not introduce
new metadata coercion, labels or invalid-input policy.

## Sparse replacement feasibility, numerical evidence and remaining gate

A scratch calculation removes only full-input `as.matrix(grna_matrix)` and sums
each selected sparse guide subgroup with `Matrix::colSums`, then uses the original
dense t(vapply()) result. The trace observes original whole-input dgC coercion at
5 by 3, followed by dense groups 2 by 3, 3 by 3 and 1 by 3. The scratch calculation
has no whole-input as.matrix event and passes sparse groups of those same shapes
to colSums. The selected Matrix methods dispatch to native `R_sparse_marginsum`;
their source bodies are recorded. This establishes a concrete sparse path,
**not measured allocation/RSS savings**.

The 200 seeded, independently copied, valid nonnegative duplicate-triplet and
metadata-rounding cases have **zero exact mismatches and zero failures at 1e-8**
on the recorded platform/dependencies. Cases include zeros, fractions, tiny
values and 2^53-scale values with metadata repetitions. Primary cases also
record exact hexadecimal rounding, for example fractional sum
`0x1.3333333333334p-1`. No universal bitwise portability claim follows from these
small observations; candidate parity on its actual platforms remains required.

Therefore **yes: the existing dense output can be retained while eliminating
full-input densification**, using the demonstrated sparse subgroup sums. The
required output still costs vectors by cells, and a selected sparse subgroup can
expand with repeated metadata. Do not replace the selector with set membership,
sum duplicate metadata only once, use a sparse output, or change supported
backends. Bounded numerical evidence exposed no remaining policy conflict.

T11's representative allocation/workspace/peak-RSS/time gate remains outstanding
before production adoption/qualification. Use fixed release builds, representative
guide/cell/vector multiplicities, cold/warm controls, a warm-up and at least three
measured repetitions; retain outliers and the required dense output footprint.
Compare median time and peak worker RSS to original under the approved gates
(+10% time, RSS increase at most max(10%,64 MiB)), investigating failures rather
than extrapolating from these tiny probes. No speed or memory improvement is
claimed here. Direct backed-vector aggregation is not added. The existing
Cell Ranger backed route delegates its optional vector map to ondisc; it is a
separate existing route, not evidence that direct R-matrix backed import has
identical aggregation or permission to extend it.

## Delivery and remaining choices

No unresolved duplicate/scientific choice remains in this bounded evidence:
B1's explicit rejection decision is recorded; B6 preserves separate coordinate
and metadata multiplicities. Production implementation and T11 resource
qualification are pending. Future parity/resource failures can still expose a
concrete issue and must be reported; this report is not approval to normalize it.

Primary machine evidence: report.json. Fixture/source hashes: fixtures.rds,
source-manifest.json, artifact-manifest.json. Scripts and fixture directories are
small, bounded, and live only under this temporary directory. The worktree and
all production sources/libraries were left unchanged.
