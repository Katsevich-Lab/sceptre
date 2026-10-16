# CL4 native closure and C1 source-review handoff

Base: `d3746f6f9aecd0dfb6bf954ce982dd81659b0026`. Exact current source hashes,
49 named retained implementation functions and their owning source contracts/
actual callers are recorded in [native-closure-coverage.json](native-closure-coverage.json).
This is a source review unit. No R process, compilation, test or resource job
was started while the independent dense/triplet resource window was active.
Runtime qualification and final integrated-package qualification are separate.

The approved [symbol dispositions](symbol-dispositions.md) govern this unit.
Public exports/formals/S4 slots, original frozen fixtures/bridge/include/configs,
stage budgets/tolerance and statistical policies remain unchanged. R import,
storage, storage-ownership tests, resource tooling/results, workflow YAML and the
root-owned dispersion/EM implementations are outside this unit.

## Coherent closure

- `run_low_level_test_full_v4` becomes `run_association_test`; four real workers,
  prepared native/index/safety tests and component record calls migrate together.
- Production `crt_index_sampler_fast` becomes `sample_crt_indices`; the unused
  slow `crt_index_sampler` is deleted rather than kept behind an alias.
- `sample_combinations_v2`/`sample_undercover_pairs_v2` become
  `sample_calibration_groups`/`sample_calibration_pairs`.
  `increment_matrix` becomes the explicit mutator
  `increment_group_indices_in_place`.
- `construct_negative_control_pairs_v2` becomes
  `construct_negative_control_pairs`, including its real calibration facade
  caller. Its metadata and low_moi locals were demonstrably unread and removed.
- `legacy_result_record` becomes `association_result_record` without changing
  its fields or arithmetic.
- Delete the unused observed/null/effect wrappers, `compute_colwise_max`,
  now-uncalled `increment_vector` and all four `check_functions.cpp` diagnostics.
- The two real tail declarations belong to `inference.h`/the inference namespace;
  the moment fit, goodness-of-fit check and outlier gate are tail-local in an
  anonymous namespace. The outlier gate now accepts a const segment it only
  reads. Delete the obsolete all-purpose shared header.
- `obtain_pointer_vector`, `load_csr_row`, loaded-vector primitives and
  `legacy_pair_eligible` retain their real responsibilities/names. No storage
  caller rename or new compatibility alias is introduced.

Both generated bindings were regenerated together using
`Rcpp::compileAttributes()` after the native owner's explicit quiet-window
release. A second generation was byte-idempotent (Rcpp 1.1.1.1.1), and all 13
changed R source/test/tool files parsed. Source/bindings contain no retired
production names or old registration. No handwritten binding patch was used.
Generated-file hashes are in the coverage JSON. Compilation, installed tests
and independent source review remain pending.

## Three active sampler comparisons, four historical configurations

The candidate map in `helper-refactor-native.R` is literally:

```r
c(fisher_yates = "fisher_yates_sampler",
  hybrid = "hybrid_fisher_iwor_sampler",
  crt_fast_production = "sample_crt_indices")
```

Each actual bank is produced by its mapped function in the **installed candidate
namespace**, then fully materialized. Expected banks come from the independent
unchanged original-source bridge. Fast CRT compares all 5,498 draws. Candidate
existence does not select coverage; missing active functions fail.

A dedicated test asserts the exact four original configurations, arguments,
source/include hashes and archived fixture hashes on every platform. The
retired slow CRT configuration has no candidate substitute: its original-source
oracle is compared with its historical full bank on the recorded platform/
toolchain. That is archive integrity, not active production parity. The original
bridge retains all four oracle dispatch cases. Historical reports are untouched.

The historical `generate-native.R` still verifies and runs only the pinned
original source. Its local original-entry runners keep that tool independent
of current renamed candidate helpers; its four configurations and numerical
recipes are unchanged. It is **not executed** by this refactor and no fixture
is regenerated. The current ownership probe covers the three retained samplers
at three controlled R pointer-boundary failures: nine cases, not twelve.

## Deliberate native contract coverage

The scoped implementation coverage is **49 named functions in 13 reviewed
native files** (ten .cpp files and three headers), plus the nine inference
interface declarations, their record/view types, and the named cis-distance
lambda. `compute_nb_size_param_functions.cpp` and `mixture_functs.cpp` are
separately root-owned/reviewed and excluded from these 49. Existing accurate
cis/count/policy contracts were inspected and retained or supplemented; no
mechanical annotation count is claimed as evidence of correctness.

| Owning file(s) | Functions covered | Responsibility and retained complexity |
| --- | ---: | --- |
| generate_samples_functions.cpp | 5 | FY ordered extraction; hybrid N+k prefix laws; CRT column-wise binomial/subset stream; fixed Boost seeds; full-bank reuse and NULL-XPtr-before-native allocation |
| low_level_full_test.cpp | 1 | Actual Rcpp entry; aligned fit arrays; one observed base conversion; borrowed const bank; observed versus synthetic cardinalities; result schema |
| inference_scores.cpp | 5 | Shape/observed-coordinate boundaries, fold ratio/SE, studentizer, selected bank segments, fixed prefixes versus entire CRT rows, degenerate arithmetic |
| inference_stages.cpp | 2 | Inclusive .02 transition, disjoint B1/B2/B3 offsets, optional SN fit and selected-segment fallback, including B3=0; real result record |
| inference_tail.cpp | 5 | Inclusive ties/add-one tails and unclamped two-sided values; moment-fit asymmetric rule; tail/outlier gates; sorting/sign reflection and selected-segment ordering |
| negative_control_functions.cpp | 7 | Guarded combinations, fixed std seed (distinct from Boost), subset rejection/row order, explicit scratch mutation, candidate stream/acceptance order and existing shared count policy |
| sparse_counts.cpp | 9 | R-rooted normalized views, one-based source/analysis domains, pure native copy/compute scope, fused stored-position counts, repeated memberships and NT-group sums |
| sparse_matrix_functs.cpp | 6 | Actual CSR row load/pointer construction, loaded-value summaries versus stored positions, zero/fractional values, first-max ties, membership multiplicity, stored-value threshold order |
| pair_construction_functs.cpp | 1 | Inclusive signed-wide distance predicate, input order, no sorted-input assumption, exact allocation without native scratch across R allocation |
| pair_eligibility.h/.cpp | 3 | Inclusive legacy statistical policy and R FALSE-dominates-NA vector boundary; this legacy label has a real policy meaning |
| sampling_safety.h | 5 | Checked integer/extent arithmetic and nonallocating final ownership transfer, without changing budgets or draw order |
| inference.h | 9 declarations/3 interface types | Shared q-by-n score/bank/result ownership and cross-translation-unit tail declarations; counted separately from implementation functions |

No new framework, statistical fallback or universal index representation is
introduced. Fixed native seeds consume no R draws; generated Rcpp RNGScope
behavior is unchanged (it can initialize a previously absent R seed), rather
than claiming that every wrapper leaves an absent `.Random.seed` absent. Stored sparse positions (including explicit zeros/fractions) are
not interchangeable with loaded-value nonzero counts (>0.5). Hybrid suffixes
can exceed a smaller current fit universe; only the selected prefix belongs to
that fit. CRT synthetic cardinalities can differ from the observed count.
These distinctions explain real retained complexity rather than obsolete
compatibility layers.

## Static evidence and pending runtime gate

The task-specific `native-closure-static.py` uses lizard 1.24.1 to enumerate
functions, then compares lexical body tokens after exactly the approved symbol
crosswalk. It confirms all 49 retained bodies match the base, signatures match
except the private read-only outlier const qualifier, exactly the ten
approved unused bodies disappear, and the three production R caller files
match after renames plus the two unread-local deletions. It also confirms the
frozen fixtures, API files, root-owned implementations and excluded import/
storage files are unchanged in this unit.

```sh
PYTHONPATH=/path/to/lizard-1.24.1 \
  python3 tools/refactor/native-closure-static.py SOURCE_ROOT NEW_OUTPUT_JSON
```

This is bounded source-body evidence. It does not prove comments correct,
namespace/linkage resolution, binding registration, ABI, runtime parity or
resource acceptance. Actual source/caller reread supplies the contract evidence;
installed tests and component/ownership probes supply the later execution gate.
Generation and source/binding absence readback are complete. Pin the coherent
source for independent senior review, then coordinate a fresh immutable release
install and focused native/full-bank/
index/calibration/context/API/saved-object checks plus the component and nine-case
ownership probes. Frozen inputs/expected files must remain byte-identical.

Root separately updates the original-versus-candidate allocation tracing map,
architecture/CI/global metrics and the historical scorecard label. Old names in
those historical selectors, original-only generator/config assertions or archived
reports are intentional evidence; they are not candidate aliases.
