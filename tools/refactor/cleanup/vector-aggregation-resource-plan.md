# B6/T11 vector aggregation acceptance protocol

Predeclared before resource execution. Candidate branch starts at 7df5ad1; the
conditional change removes full-input densification only in in-memory vector
aggregation. It retains the original dplyr metadata selector, `t(vapply())`,
dense numeric output, first-seen vector/cell order and Matrix duplicate semantics.
No new backed aggregation or output/count/rounding policy is authorized.

1. Build the candidate from an immutable Git archive using a fresh `R CMD INSTALL
   --preclean` with the same default -O2 toolchain/dependencies as original source
   21f9ea0. Preserve both libraries, source/compile logs and DLL hashes. Use source
   ondisc 1.3.5, one worker, and OMP/OPENBLAS/MKL/VECLIB thread variables all one;
   set data.table threads to one. Never replace the system ondisc installation.
2. Before large execution, replay all 28 original characterization fixtures and
   200 independently generated rounding/duplicate cases through the two installed
   namespaces in separate processes/copies. Require exact output equality and
   input nonmutation; report any tolerance comparison separately on a mismatch.
   Run independent hand-specified/public-import tests too. Do not regenerate any
   historical reference artifact or route candidate actuals through old code.
3. Generate deterministic sparse fixtures without a dense assay. Primary:
   3,200 numeric guides x 24,000 cells, 120 stored entries per guide (0.5%), four
   guides per vector -> 800 output rows. Dense input counterfactual is 614,400,000
   bytes; unchanged dense output is 153,600,000 bytes. The fixture has integral
   counts, named guides/cells and reversed metadata order. Secondary numeric
   triplet and logical CSC: 800 guides x 6,000 cells, 0.5% stored density, four
   guides/vector -> 200 dense rows (9,600,000 output bytes). Numeric triplets add
   valid duplicate coordinates, exact binary fractions, explicit zeros, reversed
   stored order and repeated metadata. Logical CSC retains explicit FALSE. Add
   actual Papalexi example guide counts when locally available, with transparently
   synthetic groups of four guides; this is not real biological vector metadata.
   All shapes, stored/nonzero densities, multiplicity, versions and input SHA256s
   are recorded. An unavailable real case is reported, never silently replaced.
4. Use one fresh child per build/case for correctness. Each child loads an
   independent sparse fixture, retains caller input, verifies nonmutation, and
   saves its complete dense result outside any performance measurement. Compare
   full result bytes/objects, dimensions, class, dimnames, values and ordering.
5. For each case run one warm-up pair and three measured fresh-process pairs,
   alternating build order original-first, candidate-first, original-first,
   candidate-first. Retain every trial. Input I/O, explicit pre-call GC and
   validation stay outside aggregation timers; automatic GC inside the helper
   remains measured. Keep input and output live at the end marker. Single-call
   timings below one second are diagnostic. For such a setting, predeclared
   timing-only batches run 20 independent loads/conversions per child (one live
   input/output at a time, explicit GC between calls outside timers), again one
   warm-up and three measured pairs. Compare the sum of the 20 conversion
   intervals, not elapsed I/O/GC. Preserve single-call RSS as the primary memory
   measurement. If even a batch is under one second, mark timing inconclusive;
   do not silently enlarge the recipe or claim a stable ratio.
6. Launch each single-process child under macOS `/usr/bin/time -l`, recording its
   whole-process OS maximum RSS in bytes. Sample the actual R PID every 20 ms
   with input-ready/aggregation-start/end markers, recording pre-call, active
   peak and post-call RSS. Pre-call and post-call markers wait for a monitor
   acknowledgement outside the conversion timer so their samples are real
   boundaries; an absent monitor stops after ten seconds. Sampling can miss brief peaks; OS whole-process RSS
   is the peak authority. No fork/shared-page aggregate applies here. Startup,
   input, required output, coercion and temporary workspace are included.
7. Run one additional allocation-profile child per build/case, separately from
   timing. Retain Rprofmem traces (64 KiB threshold), counts/total/largest traced allocations, allocations
   at least the dense-input size, stack examples and R heap high-water figures.
   This distinguishes the original full-assay allocation from the unchanged
   dense result and transpose. R allocation tracing can omit native workspace;
   RSS includes it. Active RSS minus pre-call RSS is a net increment, not an
   isolated measurement of sparse normalization/sorting workspace.
8. Before large children, record host memory and require at least 8 GiB of
   available-memory proxy (free + inactive + speculative vm_stat pages), without
   a competing heavy job. Record that this proxy is not a hard OS guarantee.
   Both builds use R_MAX_VSIZE=1536M. Stop an experiment child conservatively at
   2 GiB RSS or allocation failure, preserve all logs, and review the recipe;
   the R heap limit and sampled RSS are not hard limits on all native workspace.
   No deliberate OOM, full-input dense fixture generation or competing heavy
   work is permitted. Coordinate the whole window with the native resource owner.
9. Gates: exact structures/order/duplicate/nonmutation contracts; meaningful
   candidate median elapsed increase <=10%; candidate maximum peak RSS increase
   <=max(10% of reference peak,64 MiB). Report maximum and median per-trial OS and
   active-phase peaks, plus all traces/outliers. Investigate any crossing; retain
   targeted repeats separately and request explicit engineering exception review
   for a persistent bounded tradeoff. An unrun/unsafe/inconclusive primary gate
   is incomplete qualification, not a pass. No integration before these gates.
