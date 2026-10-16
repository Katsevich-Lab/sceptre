# Refactor resource verification

The original `benchmark.py` retains the RF0 small-screen serial diagnostic.
`resource-data.R`, `resource-benchmark.py` and `resource-worker.R` add explicit
review workloads without changing production methods or fixture expectations.

The fixed input recipe includes a compact screen; actual Papalexi2021 counts
from `sceptredata` (20,729 cells, 299 responses); deterministic synthetic sparse
counts with small NB sizes, batch-limited expression and a near-collinear
design; a guide-oriented CRT workload with repeated responses; and a backed
permutation assay with a 512 MiB dense counterfactual. The latter exceeds the
declared 256 MiB response-assay materialization budget. That budget is not an OS
memory limit, and the example does not establish that the assay cannot fit the
48 GiB Mac. Generated sparse matrices, not full dense assays, are written to
temporary input files. Real Replogle count matrices were not found in the
bounded local data inventory; the sparse difficult case is explicitly synthetic.

```sh
Rscript --vanilla tools/refactor/resource-data.R \
  tests/testthat/fixtures/refactor /private/tmp/sceptre-resource-data
R_LIBS=/private/tmp/sceptre-rf0-ondisc-library \
python3 tools/refactor/resource-benchmark.py \
  /path/to/frozen/package-library /private/tmp/sceptre-resource-data \
  tests/testthat/fixtures/refactor /private/tmp/sceptre-resource-reference \
  --label reference --source-commit VERIFIED_COMMIT
```

Each case runs cold and warm, serially and with two fork workers. Each setting
has one warm-up and three measured fresh R processes, for 80 processes per build.
Thread variables are fixed to 1. Public parameter selection retains complete
production B1/B2/B3 budgets and original sharing scopes; the JSON reports the
actual budgets and passing pairs. Preparation is outside the analysis timer.
Warm trials run discovery once before measuring the next discovery call.
`--resume` only accepts the same recorded source/input/tool/execution plan.

The process sampler reads `ps` every 0.1 seconds during the analysis phase. It
records the main R process, each observed fork descendant, maximum individual
worker RSS and synchronized aggregate tree RSS. Sampled maxima are lower bounds
and can miss short-lived peaks. Aggregate RSS counts shared pages repeatedly;
it is not proportional set size or unique physical allocation. The Mac sandbox
blocks `ps`, so execution requires the same narrow process-read escalation used
for the earlier `/usr/bin/time -l` diagnostics. No system dependency is installed.

Compare reports with `resource-compare.py`. The initial engineering review
threshold is at most 10% higher median elapsed time for reference medians of at
least 1 second, and at most the greater of 10% or 64 MiB higher sampled peak RSS.
Apply the memory comparison to main, individual worker and aggregate tree
observations, retaining their different interpretations. A threshold crossing
requires investigation and an explicit disposition; it does not justify changing
the method or accepting new scientific snapshots. Exact scientific hashes are
compared independently. Any hash difference requires numerical/identity review.

Compact timings remain diagnostic. `resource-allocations.R` separately records
R allocation counts/bytes and `glm.fit`/sampler call counts for cold/warm serial
compact work. Its Rprofmem/tracing overhead is absent from elapsed trials.
R allocation profiling excludes C++/HDF5 allocations and fork memory.

Baseline and candidate inputs, library paths, source pins, dependency/platform/
BLAS details, thread settings, scientific hashes and warning lists are retained
in JSON. Benchmark outputs are observations, not fixture regeneration. Full
transcriptome Replogle performance, actual out-of-RAM scaling and operating-system
enforced memory limits remain outside this bounded acceptance exercise.
