# Bounded ownership and fit-cache benchmark

Run from the repository root with:

```sh
Rscript tools/hardening/benchmark/benchmark-boundaries.R
```

The committed `results/` snapshot measures the current source loaded by pkgload
at frozen candidate `169c16e4ba00c687c1f03fe37e9d79a11d01cd8f`,
on local Apple arm64, R 4.6.0 and data.table 1.18.4. `metadata.R` records the UTC
run time and environment. `source-hashes.csv` and the RDS metadata identify every
R source file and this benchmark script by MD5; source stability is checked
between load and completed measurements; the recorded Git HEAD is context,
not the identity of the uncommitted implementation. An optional
`SCEPTRE_BENCH_LIBRARY` selects an installed package instead.

## Ownership boundaries

Each latency uses 20 batches of 10 serialize/unserialize ownership copies. The
allocation profile measures one additional copy, with no baseline subtraction.
Values below are the final run; CSVs retain the unrounded measurements.

| Workload | Input `object.size` bytes | Median copy, ms | Observed batch-average range, ms | Profiled allocation bytes |
| --- | ---: | ---: | ---: | ---: |
| Actual five-row result table | 4,880 | <0.1 | 0–0.1 | 11,328 |
| Actual saved record including history | 48,464 | 0.10 | 0–0.10 | 25,880 |
| Actual CRT `output_amount=3` record | 1,264,216 | 1.50 | 1.40–2.00 | 566,288 |
| Synthetic nested draw graph | 591,264 | 0.75 | 0.70–2.60 | 1,091,048 |

The synthetic graph is **not scientific output level 3**: it contains 500 rows
of 128 arbitrary normal draws (64,000 numeric values), a nested data.table
attribute, and a genuine frozen RunPlan. The actual CRT record is measured
separately. The tiny table's zero batch averages reflect clock resolution, not a
zero-cost copy.

`Rprofmem()` reports recorded R allocations, not peak resident memory, retained
memory, native allocations or whole-workflow memory. `object.size()` is an R
estimate that can count shared content repeatedly; allocation totals and object
sizes need not agree. One profile does not establish a general memory bound.
Serialization temporarily materializes a serialized buffer and an owned graph;
these measurements do not qualify arbitrarily large nested result histories.
Cleanup disables profiling/removes temporary profiles and removes fit traces on
both success and error.

## Cold and warm CRT

The real high-MOI complement-control CRT fixture has 160 source cells, 155 retained
cells and five requested pairs; its plan uses B1=499, B2=4999 and B3=0 with output
level 3. Both runs use seed 172839. Exact result values, saved fit values and
warnings agree between cold and warm runs.

| Path | Actual response-fit calls | Actual perturbation-fit calls | Three elapsed times, seconds |
| --- | ---: | ---: | --- |
| Cold, prepared screen | 4 | 2 | 0.146, 0.144, 0.154 |
| Warm, completed cached screen | 0 | 0 | 0.146, 0.154, 0.156 |

Fit-call counts are measured separately with tracing; the timings are
uninstrumented. They establish that compatible warm fits are reused. The tiny
fixture's elapsed times are noisy and effectively comparable; this is **not a
speedup claim**. The process/filesystem are already warm, there are only three
repetitions, and all public validation, draw generation, inference and ownership
work remains included. No large-screen, cross-platform, distributed-process or
peak-memory conclusion follows from this local benchmark.

Raw artifacts are `copy-boundaries.csv`, `crt-cache-times.csv`,
`crt-fit-counts.csv`, `source-hashes.csv`, `metadata.R` and
`benchmark-results.rds` under `results/`.
