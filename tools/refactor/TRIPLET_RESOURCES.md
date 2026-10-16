# Triplet conversion resource harness

This implements the approved [experiment protocol](cleanup/triplet-resource-plan.md).
It measures the existing `set_matrix_accessibility()` operation in separately
installed original/candidate packages. It does not change scientific fixtures,
production source, returned sparse semantics, or resource thresholds.

Build both immutable Git archives anew with `R CMD INSTALL --preclean` and the
same release toolchain/dependencies. Verify each with
`resource-build-provenance.py`; the original pin must be
`21f9ea098b69c444a884a49bc254de711968e3b5`. The driver rejects mismatched DLLs,
nonrelease provenance, changed tools/inputs/packages, or reused trial paths.

After the coordinating agent confirms a quiet window, run these phases in order:

```sh
python3 tools/refactor/triplet-resource.py prepare OUTPUT \
  --original-library ORIGINAL_LIBRARY --candidate-library CANDIDATE_LIBRARY \
  --original-build ORIGINAL_BUILD_JSON --candidate-build CANDIDATE_BUILD_JSON \
  --original-archive ORIGINAL_GIT_TAR --candidate-archive CANDIDATE_GIT_TAR \
  --quiet-window-confirmed
```

Use the same arguments for `correctness`, then `measure`. `summarize` reads
already completed trials and does not require a quiet window. On macOS the
execution phases require process access for `ps` and `/usr/bin/time -l`; the
restricted process sandbox cannot provide those measurements.

Preparation checks at least 8 GiB reclaimable headroom, generates the four fixed
1M/8M RDS fixtures, and compiles one small benchmark-only `CLOCK_MONOTONIC` DLL
shared by both builds. Children use the same 1.5 GiB R vector heap, four thread
environment limits of one, and one data.table thread. RSS at 2 GiB stops the
experiment without qualifying it. Small-case observations provide a conservative
scaling safety check before large conversions; extrapolation is not a measurement.

Every timing child loads its own input, performs pre-call GC, and holds input and
output through the conversion-end marker. Identical 60 ms marker padding before
and after conversion is outside the timer and gives the 20 ms RSS sampler stable
input-only and live-output observations. Five-conversion batches always reload
an independent fixture between operations, with GC and I/O outside timers.

Untimed correctness children stream input serialization hashes to files and
save sparse outputs. A separate child compares exact class, shape, dimnames and
pointer/index/value slots, without dense realization. Candidate input hashes
must remain unchanged; original mutation is recorded. Allocation-profile children
trace sorting-table creation and sorting only in their private R namespaces,
retain Rprofmem and heap high-water records, and save timestamped RSS markers.
Instrumentation is separate from primary timing/memory trials.

The results directory preserves the protocol/tool sources, build and input
hashes, every warm-up/trial/profile, full process traces, Darwin OS peak bytes,
sampled conversion/phase peaks, allocation traffic and exact output comparisons.
The JSON summary distinguishes short timing gates and every measured crossing.
Rprofmem traffic is not live workspace, and sampled RSS can miss a short peak.
No `16 bytes/nnz` estimate is treated as a measured result.

A completed measurement run is not automatically an accepted resource result.
Persistent crossings need investigation or the explicit bounded correctness-cost
review described in the protocol. This focused experiment does not replace the
final integrated sparse/backed, worker-memory, cache or consumer gates.

Sorting ownership is determined by constructors: supplied sparse slots are
borrowed and duplicated with `data.table::copy`; dense-to-triplet columns,
expanded pointer coordinates and replacement logical-one values are fresh and
can be sorted directly. This avoids arithmetic copying and refcount heuristics.
For Matrix 1.7.5 the installed dense method calls `R_matrix_as_sparse`; the
[versioned source](https://github.com/cran/Matrix/blob/1.7-5/src/coerce.c#L1495)
calls `dense_as_sparse`, which allocates triplet coordinates and values
(lines 1864, 1870, 1912). The inspected source SHA256 is
`67805fdd88ccdf06a4c1d888d55df3dec92031d27b20d5410a23cec39ab657b3`.
Independent small tests cover dense/compressed/logical caller aliases and
bitwise negative-zero/NaN-payload preservation. Allocation profiles separately
trace explicit column duplication as well as table creation/sorting.
