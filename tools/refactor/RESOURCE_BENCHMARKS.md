# Refactor resource verification

The final [acceptance report](FINAL_ACCEPTANCE.md#resources-and-explicit-dispositions)
and [raw evidence](results/final-0a9f9d8/README.md) supersede earlier checkpoint
acceptance claims below. Final coverage is the complete 160-process matrix at
2ab3f5a plus 48 targeted actual 0a9f9d8 processes and reviewed Air source carry.
Two bounded memory exceptions remain explicit; thresholds were not relaxed.

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
python3 tools/refactor/resource-build-provenance.py . \
  /path/to/fresh-git-archive /path/to/frozen/package-library \
  /path/to/clean-install.log /private/tmp/build-provenance.json \
  --commit VERIFIED_COMMIT
R_LIBS=/private/tmp/sceptre-rf0-ondisc-library \
python3 tools/refactor/resource-benchmark.py \
  /path/to/frozen/package-library /private/tmp/sceptre-resource-data \
  tests/testthat/fixtures/refactor /private/tmp/sceptre-resource-reference \
  --label reference --source-commit VERIFIED_COMMIT \
  --build-provenance /private/tmp/build-provenance.json
```

Build each package from a fresh verified Git archive using `R CMD INSTALL
--preclean`; the provenance check requires matching source blobs and a successful
release install that recompiles every native source. It rejects reused objects,
`-O0`/debug overrides, or mismatched compiler/SDK/native flags. The harness checks
the recorded native DLL hash before and after trials. This gate was added after
an earlier reference library was found to reuse development objects; that entire
batch is invalidated in `results/resource-invalidated/` and cannot establish a
speedup or satisfy acceptance.

Each case runs cold and warm, serially and with two fork workers. Each setting
has one warm-up and three measured fresh R processes, for 80 processes per build.
Thread variables are fixed to 1. Public parameter selection retains complete
production B1/B2/B3 budgets and original sharing scopes; the JSON reports the
actual budgets and passing pairs. Preparation is outside the analysis timer.
Warm trials run discovery once before measuring the next discovery call. Cold
means an empty model-fit cache; operating-system filesystem caches are not
purged. Backing files are freshly written during preparation for both builds.
`--resume` only accepts the same recorded source/input/tool/execution plan.

The process sampler pauses 0.1 seconds between `ps` snapshots during analysis;
process-listing overhead adds to the observed interval. It
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
compact work. The same separate diagnostic checks response/assignment fitting
and per-guide draw sharing for the CRT case. Its Rprofmem/tracing overhead is
absent from elapsed trials.
R allocation profiling excludes C++/HDF5 allocations and fork memory.
The profiler takes six positional arguments: library, case RDS, frozen fixture
directory, `cold` or `warm`, output JSON, and `original` or `candidate`.
The final argument selects an explicit three-sampler map for that source;
missing functions fail tracing rather than silently reducing coverage. The
archived original slow CRT sampler is not an active production producer.

Baseline and candidate inputs, library paths, source pins, dependency/platform/
BLAS details, thread settings, scientific hashes and warning lists are retained
in JSON. `resource-environment.R` also records resolved dependency versions,
paths and compiled-library hashes for each frozen library. Benchmark outputs
are observations, not fixture regeneration. Full
transcriptome Replogle performance, actual out-of-RAM scaling and operating-system
enforced memory limits remain outside this bounded acceptance exercise.

## Recorded release comparison

The bounded result is in `results/resource-release/acceptance.json`, with the
complete batches, clean-build source/flag/DLL proof, dependency sidecars, input
recipes, allocation/call counts and timestamped repeat reports beside it.
Reference source is `21f9ea0`; candidate production source is `28df390`.

The primary batch ran 160 fresh processes: 120 measured repetitions and 40
warm-ups across five cases and 20 settings per build. All 80 corresponding trial
pairs match exactly for scientific hashes, budgets, passing pairs, input hashes,
warnings and runtime configuration. Compiler/SDK/flags and 13 dependency packages
match. All 11 settings with reference medians of at least one second pass the
10% elapsed review threshold; the largest candidate/reference median ratio is
1.017. Shorter timings remain diagnostic.

Main and observed individual-worker peaks pass the specified memory budget.
One initial Papalexi cold/two-worker aggregate snapshot crossed its additional
conservative aggregate threshold: candidate 1,170 MiB versus reference maximum
1,039 MiB. It remains visible in the primary comparison. A fresh timestamped
paired repeat (warm-up plus three measured trials each) observed maxima of
1,044 MiB versus 1,030 MiB, within budget, with exactly two fork workers. Peak
snapshots occurred about 1.33–1.37 seconds after the first analysis snapshot.
No persistent regression was identified; the precise original cause remains
unknown. R heap/COW/sampling variability is a possible explanation, not a
measured cause. Thresholds and production code were unchanged. An additional
four-process reference probe is also retained separately.

The separate 12-process allocation/call-count diagnostic matches scientific
hashes and fit/bank counts exactly. Compact cold/warm work performs three/zero
`glm.fit` calls and one full-budget permutation bank. CRT cold/warm work performs
72/eight `glm.fit` calls and eight guide-specific banks, preserving response-fit
and draw sharing. Fast warm sparse/backed work performs zero fits and one
permutation bank. Compact fork peaks can finish between RSS samples; their
missing observations are explicit rather than recorded as zero memory.

These results support the selected bounded workloads. They do not establish
full Replogle-scale or out-of-RAM performance, and aggregate RSS does not measure
unique physical memory. No speedup is claimed.
