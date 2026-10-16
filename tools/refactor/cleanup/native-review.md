# Native inventory input to the cleanup plan

Source: c7d62502757f2a17741e2ef18b9e4d006bf063dc. Independent senior
engineer inspected all 14 handwritten implementation files and three headers,
their R callers and the bounded actual-consumer inventory. No code was edited,
tests run or unsafe examples executed. This is an input audit; formal review
of CLEANUP_PLAN.md is a separate checkpoint.

| Source owner | Actual responsibility/caller | Register disposition |
| --- | --- | --- |
| generate_samples_functions.cpp | Permutation/CRT draw banks; preprocessing and association workers | S1 unsafe M greater than universe; S2 integer/range/bank preconditions; characterize the slow CRT sampler's test-only status separately from production kernels |
| compute_nb_size_param_functions.cpp | Dispersion estimate consumed by fitting.R | B5: numerical fallback arithmetic is real behavior; catch-all can mask nonnumerical failures and method labeling is premature |
| mixture_functs.cpp | Native EM consumed by mixture_model_functs.R | Preserve active model, initialization and convergence/fallback behavior; remove no algorithm solely on appearance |
| negative_control_functions.cpp | Calibration combinations/candidate acceptance plus a misplaced group-count kernel | S2 cardinalities/zero-rate narrowing; I2 mutation; move count ownership out of candidate sampling without changing candidate stream |
| pair_construction_functs.cpp | Public cis pair constructor | B2: unsorted positions can omit valid responses; interval arithmetic/NA contract also needs an explicit boundary |
| pair_eligibility.h/.cpp | One legacy eligibility predicate used by vector QC/calibration | Retain shared owner, complete/NA semantics and corresponding caller preconditions |
| pairwise_qc.cpp | Fused sparse pair and NT summaries via storage | I1: eliminate unnecessary NT-local round-trip, retain stored-entry and NT-group-sum semantics/order |
| shared_low_level_functions.h/.cpp | Stored-position loader plus unrelated tail declarations | D2: move loader into count ownership and tail declarations into inference; no all-purpose shared module |
| inference.h, inference_scores.cpp, inference_stages.cpp, inference_tail.cpp, low_level_full_test.cpp | Actual score/effect/tail stages and Rcpp record boundary | D1/D2: delete three unused global wrappers, retain actual R entry boundary; S2: skew-normal check requires at least 100 draws (maximum p=.99); ordinary 4,999-draw stage is safe |
| sparse_matrix_functs.cpp | Cell summaries, guide counts, thresholds and scratch offsets | I2: cell covariates decrement mt_gene_idxs; D1: unexported compute_colwise_max is unused; scratch-offset primitive has actual importer/getter callers, so do not remove blindly |
| check_functions.cpp | Unexported diagnostic functions and an unused test() | D1: no runtime/consumer caller, remove from production |

Material evidence from this audit:

- S1: preprocessing_functions.R chooses M by summing NT sizes. Two overlapping
  groups of size eight in a ten-cell universe give M=16. draw_wor_sample then
  indexes below zero. This was derived from source, not executed.
- S2: missing calibration-group bounds can reach unchecked combinations/draw
  loops; n_groups*n_genes is signed int before allocation; a zero p_hat can be
  narrowed from a nonfinite group-count estimate. The tail minimum is 100, not
  200: the current loop's largest probability is .99.
- B1: a coordinate table with cell_idx 0,1,0 reaches a pointer vector ending at
  two with three stored values; imported grouping assumptions are not established.
- B2: positions 90,1000,110, midpoint 100 and distance 20 omit the third response
  because native construction exits after its first interval block.
- B5: estimate_theta catches every C++ exception over MLE/MM. Allocation/system
  errors can become apparently usable estimates; method is set before fallback,
  and the R caller discards it. Changing exception policy needs a specific
  characterization/disposition, not an automatic refactor claim.

Unused global wrappers verified in production/tests/tool source:
compute_observed_full_statistic_v2, compute_null_full_statistics,
estimate_log_fold_change_v2. No package/pinned-consumer caller was found.
The actual Rcpp full-test entry remains a real boundary until its callers are
migrated coherently. Public NAMESPACE identity is independent of private names.

Per-file hashes and parser-derived function lists are in review-coverage.json.
Static reachability is bounded; arbitrary external use of undocumented symbols
was not surveyed. All safety/behavior items require the formal plan disposition.
