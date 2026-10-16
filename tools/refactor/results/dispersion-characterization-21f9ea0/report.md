# Dispersion numerical-exception characterization

Status: bounded read-only characterization; no production change or universal
exception/portability claim. Source is original 21f9ea098b69c444a884a49bc254de711968e3b5.

In R 4.6.0 aarch64 with BH 1.90.0.1 and Rcpp 1.1.1.1.1, 107 independently recorded
aligned-input cases returned exactly the original installed estimate/method
record: 11 completed without an exception, 61 caught std::domain_error and 35
caught std::overflow_error. The run includes seven hand-specified cases and a
seeded 100-case bounded extreme-range sweep. Inputs and raw nonfinite results
are retained in report.rds; report.json and run.R provide readable reproduction.

The instrumented source changes only the entry name/export annotation and the
empty catch's exception-type/message capture. Mathematical calls, loop/order,
returned state and fallback decisions are unchanged. It is compiled separately
from the immutable installed original; candidate code is not an expectation
source. Actual digamma/trigamma calls at zero, tiny positive, ordinary, infinity
and NaN arguments independently demonstrate domain and overflow exception paths.
All numerical exceptions observed in the complete estimator occur before the
MLE result is committed, preserving the pilot estimate/method 3. This is an
observation on these cases, not a claim that every imaginable failure was sampled.

Executed headers explain the distinction. Digamma/trigamma pole handling raises
std::domain_error; their checked narrowing/overflow raises std::overflow_error.
BH's generic policy also defines evaluation/rounding exception classes, but no
such path was observed here; fixed-precision special-function dispatch and
platform-specific paths must be checked before final catch adoption. MM uses
existing ordinary/Rcpp arithmetic rather than these special functions. A plain
catch(std::exception) would still mask unrelated allocation/system failures.

Proposed bounded correction, subject to independent review: catch the concrete
numerical exception categories required by the actual supported BH dispatch,
leaving their current estimate/method state and MLE/MM/pilot arithmetic unchanged;
propagate unrelated allocation/system errors. Controlled test-only fault injection
at MLE and MM call boundaries should verify propagation without real OOM or new
production hooks. Preserve any existing partial-state numerical path established
by review; select no new estimator/fallback/method labeling. Final-source and
cross-platform regressions remain necessary. The study is evidence for T12,
not qualification of a replacement catch policy.

The initial sourceCpp export-discovery issue is preserved in
failed-initial-loader.log. It ran no estimator assertions; suppression of the
included entry's export annotation fixed the test wrapper. The final 107-case
run completed successfully with all original comparisons exact.
