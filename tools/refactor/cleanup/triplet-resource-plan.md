# Triplet ownership correction: resource experiment plan

Status: **planning only; not executed and not qualified**. This closes the O1/C02
measurement-plan gap identified by the independent cleanup-plan review. The
underlying source audit stays pinned to
`c7d62502757f2a17741e2ef18b9e4d006bf063dc`; that pin still has the mutation. Run this
protocol only after the overall plan checkpoint permits the ownership change.
No unsafe native-boundary examples are part of this experiment.

## Comparison and immutable builds

Compare the original sceptre source
`21f9ea098b69c444a884a49bc254de711968e3b5` with the accepted candidate commit that
contains the ownership correction. Record the candidate's full commit before
execution. Build **both anew** from verified clean Git archives into separate
immutable temporary libraries with `R CMD INSTALL --preclean`, matched release
`-O2` compiler configuration and complete install logs. Do not reuse pkgload/debug
objects or overwrite system packages. Old performance reports/libraries are
historical evidence, not qualification of this new path.

The manifest records source/archive hashes, exact install commands, effective
compiler flags, package/DLL hashes before and after, R/platform, dependency
versions and paths, Matrix/data.table DLL hashes, BLAS, CPU and memory. Both
children use the same R and dependency libraries, one R process and no fork
workers. Set and verify `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`,
`MKL_NUM_THREADS=1`, `VECLIB_MAXIMUM_THREADS=1`, and
`data.table::setDTthreads(1)`. Keep the machine free of other heavy jobs. The
conversion calls are the same existing private `set_matrix_accessibility()`
entry point in the respective installed builds; do not compare different
constructors or alter the mathematical operation.

## Fixed inputs, generated outside measured children

Use these eight conversion cases in each build:

| Shape | Stored entries | Representation | Output orientation |
| --- | ---: | --- | --- |
| 20,000 features by 50,000 cells | 1,000,000 | numeric `dgTMatrix`, logical `lgTMatrix` | row and column, separately |
| 20,000 features by 50,000 cells | 8,000,000 | numeric `dgTMatrix`, logical `lgTMatrix` | row and column, separately |

Density is 0.1% and 0.8%. An eight-byte dense assay would require 8,000,000,000
bytes; **never allocate it**. The smaller case establishes scaling, and the larger
case makes owned copies and sorting workspace meaningful without a dense input.
These are synthetic conversion workloads, not claims about fitting a biological
screen. The existing general sparse/backed resource suite remains a separate gate.

Generate one immutable RDS fixture for each size/representation, using the same
coordinate recipe for numeric and logical forms:

1. Set `RNGkind("Mersenne-Twister", "Inversion", "Rejection")` and seed `271828`
   separately for each size. Let `k` be doubles `1:nnz`, and
   `q = (104729 * k + 37) %% 1000000000`. The multiplier is coprime to the grid
   size, so these positions are unique at the selected cardinalities; all
   arithmetic is below the exact-integer limit of binary doubles.
2. Let `perm = sample.int(nnz, nnz, replace = FALSE)`. Set zero-based
   `i = as.integer(q[perm] %% 20000)` and
   `j = as.integer(q[perm] %/% 20000)`. These are valid, deliberately unsorted
   triplets without duplicate coordinates. Validate uniqueness/ranges once in
   fixture preparation; do not create a full-grid index/permutation.
3. Numeric `x = ((k[perm] %% 97) + 1) / 7`, with
   `x[k[perm] %% 19 == 0] = 0`. This includes explicit stored zeros and fractions.
   Logical `x = (k[perm] %% 7 != 0)` includes explicit stored `FALSE`. Preserve the
   characterized conversion of every stored logical entry to numeric one; this
   experiment does not redefine that semantic rule.
4. Use feature names `r00001` through `r20000` and cell names `c00001` through
   `c50000`. Construct the appropriate S4 triplet directly, validate it and save
   with compression disabled. Record shape, stored count, zero/FALSE count,
   range, recipe/RNG/R version and hashes of fixtures and individual slots.

The fixtures belong to this resource recipe. Never regenerate or overwrite RF0
scientific goldens. Duplicated-coordinate policy is independently characterized;
this experiment does not coalesce or redefine it.

## Independent correctness and ownership checks

Run separate, untimed correctness children for both builds and all eight cases.
Each child reads its own fresh copy of the immutable fixture; original and
candidate never share an in-memory object. Do not feed the candidate a matrix
already reordered by the original. Capture streaming serialization hashes of the
input before and after conversion, and save the returned sparse result.

A separate comparison process checks exact class, dimensions, both dimname
vectors, pointer/index/value slots, stored-entry count, and ordering between
original and candidate returned results. No `as.matrix()` or dense comparison is
allowed. Check numeric zeros/fractions and stored logical `FALSE` explicitly.
Compare against **original sceptre**, whose column-name behavior is the relevant
contract; do not substitute ondisc's private conversion, which drops column names.

Candidate caller-input serialization must remain identical. Original input
mutation is recorded as the reproducer, not required to pass the new nonmutation
contract. The returned outputs must still be exactly equal. Independent small
regressions also cover both orientations; these large checks add resource-scale
ownership evidence without making legacy mutation an expected candidate result.

## Paired timing and total/phase memory measurement

For each size/representation/orientation setting:

- Run one warm-up child per build, discard it from summaries but retain its log.
  Then run **four measured fresh-process pairs**, with build order original-first,
  candidate-first, candidate-first, original-first. Run one child at a time.
  This is 16 warm-ups and 64 measured children over the eight settings.
- Each child loads one independent input copy before the conversion timer, fixes
  threads, performs the same explicit pre-call `gc()`, and keeps the input and
  returned output live through the end-of-conversion marker. Time one conversion
  with monotonic elapsed and user/system CPU. Exclude input I/O and later
  correctness serialization from the conversion timer. Performance children use
  inexpensive dimension/stored-count checks; large serialization/output equality
  is performed in the separate correctness runs above.
- Launch each child through `/usr/bin/time -l` (or the platform's documented
  equivalent) to obtain the **maximum RSS of the complete child**, including
  startup, input loading, conversion, owned output and sorting workspace. On
  macOS record the reported byte unit explicitly; do not mix KiB and bytes.
- Sample the child's RSS every 20 ms with timestamped load-start,
  input-ready/pre-GC, conversion-start, conversion-end and exit markers. Report
  pre-conversion RSS, the sampled conversion-window peak, post-conversion RSS,
  and OS whole-process maximum. A sampled trace can miss a short peak; the OS
  maximum remains the total-peak authority. There are no fork/shared-page sums
  in this single-process experiment.
- Profile one additional instrumented child per build/setting, separately from
  timing trials. Retain R allocation traces and R heap high-water statistics;
  mark entry/exit of sorting-table creation and `data.table::setorderv()` in
  this profiling child, without adding instrumentation to production source.
  Correlate these phases with sampled RSS. Report the observed total conversion
  workspace, including ownership copies and sorting temporaries. R allocation
  traces may omit native allocations; the RSS trace/OS peak includes them.
  `conversion_peak - pre_conversion_RSS` is a net RSS increment, not an exact
  isolated sorting allocator measurement. Report that limit rather than claiming
  that 16 bytes/nnz accounts for all workspace.

Report all individual trials, medians/ranges, whole-process and conversion-phase
peaks, and reference/candidate deltas. The 1M case may be too short for an elapsed
ratio. Predeclare a one-second median conversion threshold for a meaningful
single-call timing comparison. If even the 8M case is shorter, add four paired
**timing-only batch** trials: five conversions per child, each using a newly
loaded independent fixture, one live input/output at a time, with fixed `gc()`
between conversions outside timers. Sum only the five conversion intervals.
Do not repeatedly convert the original's already sorted object. Report batch
results separately; primary fresh single-conversion memory measurements remain
unchanged. If timing is still too short, label that gate inconclusive and use
operation/allocation evidence; do not claim an unstable timing pass or silently
enlarge the workload.

## Practical memory bound and stop conditions

Before large trials, record host available memory and ensure at least 8 GiB of
headroom with no competing heavy job. Set the same 1.5 GiB R vector-heap limit
in both builds (`R_MAX_VSIZE=1536M`) and monitor one process at a time. Use a
conservative RSS safety stop at 2 GiB; terminate the current experiment child
and preserve logs if it reaches that threshold or reports allocation failure.
The R vector limit does not cap data.table/native workspace, and an RSS sampler
is not a hard OS memory limit. The headroom check and conservative stop are
therefore both required. Run the 1M cases first; do not launch 8M cases if their
observed scaling threatens the available host budget. No OOM probe is permitted.

A safety stop, missing headroom or unrun large case is an **incomplete resource
qualification**, not a pass or an automatic exemption. Investigate or review a
revised bounded recipe before continuing; record the deviation rather than
quietly shrinking the input until the result looks acceptable.

## Gates, exception review and durable evidence

The existing general gates remain: for meaningful timings, candidate median
elapsed increase at most 10%; candidate peak RSS increase at most the greater of
10% of reference peak and 64 MiB. Show both the maximum and median per-trial RSS
peaks; assess whole-process maximums and conversion-phase peaks separately, and
investigate a persistent crossing in either view. Preserve every initial
crossing/outlier. A repeated paired run may assess persistence, but must remain
separate and must not replace the original record.

If the ownership correction crosses a general gate, choose and record one of:

1. Investigate and reduce avoidable copies/sorting workspace, keeping candidate
   nonmutation and exact returned semantics, then remeasure the changed immutable
   build under the same recipe; or
2. Request an **explicit reviewed bounded correctness-cost exception**. Supply
   per-case elapsed/RSS/workspace deltas and bytes-per-nnz scaling, the demonstrated
   caller-mutation fix, alternatives considered, practical memory limit, remaining
   risks and reviewer disposition. This is an accepted exception only after that
   review; it is not a general-gate pass and not an automatic allowance of
   approximately 16 bytes/nnz.

No threshold is silently relaxed. No exception to this focused triplet case
waives the ordinary end-to-end sparse/backed, worker-memory or fit/bank sharing
resource gates. Earlier release reports remain historical and do not qualify the
new sorting path. Do not promise a speedup.

Persist a small JSON/Markdown report containing protocol/build/input hashes,
thread/dependency readback, process/phase units, every trial/warm-up/profile,
input mutation and returned-output equality results, measured summaries,
threshold comparisons, safety/incomplete outcomes and any explicit exception
review. Keep large fixture/output/profile files in a bounded temporary results
directory with hashes and recipe; do not commit them as scientific fixtures.
