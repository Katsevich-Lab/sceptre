# Native inference components

`src/low_level_full_test.cpp` keeps the existing Rcpp entrypoint and native
helper signatures. It unwraps the existing worker-local draw pointer and calls
ordinary private components declared in `src/inference.h`:

| Owner | Responsibility |
| --- | --- |
| `inference_scores.cpp` | Effect estimate, observed score, and a requested segment of null scores. Observed treatment positions are one-based; draw-bank positions are already zero-based. |
| `inference_tail.cpp` | Existing empirical and skew-normal calculations, including fit checks and in-place ordering of an owned null-score vector. |
| `inference_stages.cpp` | Existing stage threshold, segment offsets, fallback decisions, and conversion to the legacy return record. |
| `shared_low_level_functions.cpp` | Sparse nonzero loading used by pairwise QC, separate from inference. |

`TailResult` owns only the selected null-score segment. The input draw bank is a
const reference with the caller's worker-local lifetime. Synthetic banks, fitted
response arrays, and native pointers are not cached or serialized by these
components. There is no common method registry or preparation requirement for a
future inference procedure.

Run the direct component checks from a source checkout:

```sh
Rscript tools/refactor/test-native-components.R
```

This developer-only bridge compiles the actual candidate modules and exercises
their independent inputs. Its expected staged outputs come from the immutable
RF0 fixtures. The installed-package `refactor-native` tests independently check
the compatibility entrypoint and the frozen same-environment sampler oracle.
Neither command regenerates fixture expectations.
