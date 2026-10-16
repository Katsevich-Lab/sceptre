# B6/T11 sparse vector aggregation evidence

Measured source/tool candidate: `88ed37cd199c642edd0ba0c008ab507bab01d10f`;
base: `7df5ad1b164b00ae45f91a0ded3085c45764a8ca`.
Original release: `21f9ea098b69c444a884a49bc254de711968e3b5`.
Integrate the ordered source/tool commits `3a013c6` then `88ed37c` together.
This evidence concerns B6; the final integrated package still requires its own
qualification and independent source/resource review.

The helper converts a sparse guide matrix to Matrix `RsparseMatrix` once,
then keeps the original dplyr metadata selector and `t(vapply())` dense numeric
assembly. Nonsparse inputs retain the original `as.matrix()` boundary. Numeric
and logical Matrix semantics, valid triplet duplicate aggregation, repeated
metadata selections, first-seen vector order, source-cell order, dimnames,
rounding and input nonmutation are preserved. Storage's logical stored-position
conversion is deliberately not involved. There is no new backed aggregation.

## Exact and independent checks

Two separate installed namespaces produced exactly identical complete outputs
for all 28 immutable characterization fixtures and 200 independently generated
valid triplet/rounding/duplicate cases. Both inputs remained unchanged. The
frozen fixture SHA256 is
`f53c93eabbb30a53d5fed01788ee8ef0ef4447b94ed308540410e21baa624e2d`;
it was not regenerated. The two replay RDS files are observed experiment
results, not golden test expectations. Independent hand-specified helper and
public-import tests, existing coordinate imports and legacy imports passed
**15 tests / 134 expectations**, with no failure, error, warning or skip.

The large-count/fraction test explicitly retains observed original double
accumulation: increments of one following `2^53` are lost, and `.1+.2+.3`
rounds above `.6`. No new count convention or compensation is introduced.

All four resource inputs also produced byte-identical complete dense results
in independent original/candidate correctness processes. All 80 children have
matching dependency versions, paths, DESCRIPTION and DLL hashes, R/BLAS/RNG
settings, one-worker/thread settings, file/object input hashes and within-case
output structures. Every input is unchanged and every warning list is empty.
Complete outputs are serialized only in the correctness children; timing
children verify structure and input nonmutation outside the conversion timer.

## Predeclared bounded experiment

The [protocol](../../cleanup/vector-aggregation-resource-plan.md) was committed
before execution. Both immutable archives were freshly compiled with identical
release `-O2` prefixes; build/source/compiler/SDK/DLL/log pins are retained.
R 4.6.0, Matrix 1.7.5 and source ondisc 1.3.5 were used; the system ondisc was
not replaced. All four thread variables and data.table threads were one;
there was one R worker and no forks. Root and the native worker explicitly
held competing R jobs throughout this measurement window; unrelated host
activity cannot be excluded.

Each setting had one warmup plus three measured fresh-process pairs, with
alternating build order. Subsecond single calls triggered the predeclared
20-conversion timing-only batches, again one warmup plus three measured pairs.
Input loading, explicit GC, monitor acknowledgements and validation were outside
conversion timers. All 80 children are retained: 8 complete-output checks,
32 single-call children, 32 batch children and 8 separate allocation profiles.
There was one complete resource attempt, with no failed/dropped trial, safety
stop or threshold crossing.

No full dense fixture was generated. Primary shape was 3,200 guides x 24,000
cells, 384,000 stored counts (0.5%), grouped four per vector into an unchanged
800 x 24,000 dense result. The former dense input payload is 614,400,000 bytes;
the required output payload is 153,600,000 bytes. Secondary shapes were
800 x 6,000 numeric triplets (valid duplicate coordinates, fractions, stored
zeros and repeated metadata) and logical CSC (stored FALSE). Actual Papalexi
guide counts were 110 x 20,729, with **synthetic** groups of four (last group two)
into 28 vectors; these are not biological vector annotations. The numeric
triplet has 24,240 stored entries after duplicates, i.e. density 0.00505;
input-manifest JSON prints density rounded, while counts and report.json retain
the exact fraction.

Each child required at least 8 GiB free+inactive+speculative vm_stat proxy
headroom. Both builds used `R_MAX_VSIZE=1536M`; the monitor stopped only its own
child at 2 GiB sampled RSS. The proxy and sampled stop are not hard limits on
native allocation. All children completed within those declared bounds.

| Setting | Original median 20 calls (s) | Candidate (s) | Timing gate | Original max measured single RSS (MiB) | Candidate (MiB) | RSS gate |
| --- | ---: | ---: | --- | ---: | ---: | --- |
| Large numeric CSC | 20.912 | 12.353 | pass | 1374.14 | 613.61 | pass |
| Numeric triplet + metadata repeats | 2.043 | 1.686 | pass | 326.45 | 279.42 | pass |
| Logical CSC + stored FALSE | 1.982 | 1.719 | pass | 296.09 | 280.48 | pass |
| Actual Papalexi guides, synthetic vectors | 0.673 | 0.467 | inconclusive: under 1 s | 293.94 | 271.72 | pass |

The gates were a meaningful median elapsed increase <=10% and a peak RSS
increase <=max(10%,64 MiB). Actual-data timing remains explicitly inconclusive;
there is no actual-data speed claim or enlarged post-hoc batch. The primary
and both meaningful secondary timing gates pass. Every OS and sampled active
RSS gate passes. Warmups and all individual measured values remain in the raw
evidence, including warmup peaks higher than some measured candidate peaks.

## Allocation and memory interpretation

In the primary profile, the original performs a **614,400,048-byte** allocation
through `as.matrix()`. The candidate has no allocation at least that input size.
Both retain the required **153,600,048-byte** dense output assembly and a second
same-size transpose. The largest remaining candidate allocation besides those
is 3,072,048 bytes at sparse conversion. For allocations >=64 KiB, primary
original/candidate totals are 3,208/809 allocations and
1,921,036,704/466,329,552 allocated bytes. These are cumulative allocation
volumes, not live RSS or a claim about bytes per stored entry.

All four original profiles contain one full-input-sized allocation; every
candidate profile contains none. Complete Rprofmem traces and stack examples
are retained. Rprofmem omits smaller allocations and may omit native malloc;
OS RSS includes process startup, live inputs/results and temporary workspace.
Twenty-millisecond actual-R-PID samples can miss brief peaks, so macOS
`time -l` whole-child maximum RSS in bytes is the primary memory authority.
Active-minus-pre-call RSS is only a net increment, not isolated normalization
workspace. No fork/shared-page aggregate applies.

## Evidence and reproduction

- `report.json`: source-pinned conclusions, all case medians/maxima, allocation
  summaries, exact output hashes and limitations.
- `all-trials.json.gz`: all 80 full child records, phase timestamps, individual
  conversion intervals, OS time text, headroom, full sampled RSS and environment.
- `allocation-traces.json.gz`: all eight complete allocation traces.
- `attempt-file-manifest.json`: SHA256/size of every file from the full attempt,
  including large dense outputs retained under
  `/private/tmp/sceptre-b6-resource-88ed37c-run1` rather than committed.
- `candidate-build.json`, `reference-build.json`, `input-manifest.json`: release
  and deterministic fixture provenance; original and candidate DLLs were
  verified unchanged after execution.
- `pre-resource-report.json`, `focused-runner.R`, `focused.log`,
  `replay-reference.rds`, `replay-candidate.rds`: independent installed exact
  qualification. The replay runner hashes and input recipe are in the tools.
- `harness-smoke-attempts.json.gz`: preserved first sandbox-only smoke (ps was
  blocked) and passing narrowly escalated original-only repeat; no resource
  acceptance is inferred from the tiny smoke.

From the candidate source pin, create fresh archive/library/output paths rather
than overwriting existing artifacts. Rebuild with `R CMD INSTALL --preclean`,
verify with `resource-build-provenance.py`, and replay each installed namespace
with `aggregation-characterization.R`. Compare the resulting `outputs` and
`generated_input_sha256` fields exactly. Run the focused test files in
`focused-runner.R` against the installed candidate. Under a coordinated quiet
window, use `aggregation-data.R NEW_DATA_DIRECTORY`, then:

```sh
R_LIBS=/private/tmp/sceptre-rf0-ondisc-library \
  python3 tools/refactor/aggregation-benchmark.py \
  REFERENCE_LIBRARY CANDIDATE_LIBRARY REFERENCE_BUILD_JSON CANDIDATE_BUILD_JSON \
  NEW_DATA_DIRECTORY NEW_OUTPUT_DIRECTORY
```

The runner pins all four thread variables and R heap limit itself and requires
macOS ps/time observation permission. It rejects an existing output directory
and verifies input and installed-DLL hashes. Preserve any failed attempt and
investigate before a separately reviewed repeat. No baseline, golden output,
public API, scientific method, backed capability or duplicate policy was changed.
