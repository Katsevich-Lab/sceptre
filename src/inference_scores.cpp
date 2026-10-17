#include "inference.h"
#include <cmath>
#include "sampling_safety.h"

using namespace Rcpp;

namespace sceptre {
namespace inference {

// Check that a and w have one entry per cell and D has one column per cell.
// The R entry separately checks counts, fitted means and treatment-group size.
// D may have zero rows. This check does not assess design rank or resampled
// cell positions; those are the fitter's and sampler's responsibilities.
void check_score_dimensions(const NumericVector& a, const NumericVector& w,
                            const NumericMatrix& D) {
  if (a.size() != w.size() || a.size() != D.ncol())
    stop("Score arrays a/w and D columns must have aligned cell dimensions.");
}

// Copy observed treatment positions and convert from R numbering (1, 2, ...)
// to C++ numbering (0, 1, ...), keeping order and repeats. Reject missing or
// nonpositive positions here; the R entry checks the upper limit. Never subtract
// from the caller's original vector, which could change saved assignments.
IntegerVector observed_positions_zero_based(const IntegerVector& positions) {
  IntegerVector zero_based(positions.size());
  for (R_xlen_t i = 0; i < positions.size(); ++i) {
    if (positions[i] < 1) stop("Observed treatment positions must be positive and not NA.");
    zero_based[i] = positions[i] - 1;
  }
  return zero_based;
}

// Estimate fold change as the treated cells' total counts divided by their total
// fitted means, with the current residual-based standard error. Both are on the
// fold-change scale, not the log scale. y and mu use the same selected cells;
// trt_idxs gives their positions starting at 0. Repeated positions count again.
// Keep summation order and ordinary NaN/Inf results for zero or undefined ratios.
EffectEstimate estimate_effect(const NumericVector& y, const NumericVector& mu, const IntegerVector& trt_idxs, int n_trt) {
  double sum_y = 0, sum_mu = 0, sum_top = 0, curr_top, n_trt_double = (double) n_trt;
  for (int i = 0; i < n_trt; i ++) {
    sum_y += y[trt_idxs[i]];
    sum_mu += mu[trt_idxs[i]];
  }
  double fc = (sum_y/n_trt_double)/(sum_mu/n_trt_double);
  for (int i = 0; i < n_trt; i ++) {
    curr_top = y[trt_idxs[i]] - fc * mu[trt_idxs[i]];
    sum_top += curr_top * curr_top;
  }
  double se = sqrt(sum_top/(sum_mu * sum_mu));
  return EffectEstimate {fc, se};
}


// Compute the standardized score for treated cells T:
//   sum(a[T]) / sqrt(sum(w[T]) - ||sum(D[, T])||^2).
// D has cells in columns and eigenvectors in rows, in R's eigen() order. Its row
// count equals the number of design coefficients. T starts at 0 and has already
// been checked. Keep repeated positions and summation order. A zero or negative
// variance keeps the existing NaN/Inf behavior; do not add a numerical fallback.
double observed_score(const NumericVector& a, const NumericVector& w, const NumericMatrix& D, const IntegerVector& trt_idxs) {
  double lower_right = 0, lower_left = 0, top = 0, inner_sum;
  int D_nrow = D.nrow();

  // Sum each covariate-adjustment row over the treated cells.
  for (int i = 0; i < D_nrow; i ++) {
    inner_sum = 0;
    for (int j = 0; j < trt_idxs.size(); j ++) {
      inner_sum += D(i, trt_idxs[j]);
    }
    lower_right += inner_sum * inner_sum;
  }

  // Sum the treatment weights (variance term) and residuals (numerator).
  for (int j = 0; j < trt_idxs.size(); j ++) {
    top += a[trt_idxs[j]];
    lower_left += w[trt_idxs[j]];
  }

  // finally, compute the z-score
  return(top/sqrt(lower_left - lower_right));
}


// Calculate B null scores using draws start_pos through start_pos+B-1.
// Each draw lists cell positions starting at 0. For permutations, read only the
// first n_trt positions. For CRT (use_all_cells), use the entire draw, whose size
// may differ from the observed treatment group or be zero. The sampler makes
// positions valid; check the selected draw lengths and requested range here.
// Hybrid draws can contain unused positions belonging to a larger fit, so do
// not check or read the unused remainder. Keep the full draws unchanged for
// other genes and stages. Check array dimensions at the R entry, before creating
// C++ result buffers: even D.ncol() can use Rcpp protection that raises an R error.
std::vector<double> null_score_segment(const NumericVector& a, const NumericVector& w, const NumericMatrix& D, int start_pos, int B, int n_trt, bool use_all_cells, const DrawBank& synthetic_idxs) {
  sceptre::sampling::nonnegative(n_trt, "n_trt");
  const int end_pos = sceptre::sampling::sum(start_pos, B, "Draw segment end");
  if (static_cast<size_t>(end_pos) > synthetic_idxs.size())
    stop("Selected draw segment exceeds the bank length.");
  // initialize variables
  const int* curr_vect;
  // Use Rcpp's already cached data pointer; this call does not allocate.
  const double* D_data = D.begin();
  std::vector<double> out(B);
  double lower_right = 0, lower_left = 0, top = 0, inner_sum;
  int D_nrow = D.nrow(), idx;

  // iterate over the B tests
  for (int k = start_pos; k < end_pos; k ++) {
    if (synthetic_idxs[k].size() > static_cast<size_t>(std::numeric_limits<int>::max()))
      stop("Draw row exceeds native integer capacity.");
    if (use_all_cells) n_trt = synthetic_idxs[k].size();
    else if (synthetic_idxs[k].size() < static_cast<size_t>(n_trt))
      stop("Selected draw row is shorter than the requested treatment prefix.");
    curr_vect = synthetic_idxs[k].data();
    lower_right = 0;
    inner_sum = 0;

    // Sum each covariate-adjustment row over the treated cells.
    for (int i = 0; i < D_nrow; i ++) {
      inner_sum = 0;
      for (int j = 0; j < n_trt; j ++) {
        inner_sum += D_data[i + static_cast<R_xlen_t>(D_nrow) * curr_vect[j]];
      }
      lower_right += inner_sum * inner_sum;
    }

    // Sum the treatment weights (variance term) and residuals (numerator).
    lower_left = 0;
    top = 0;
    for (int j = 0; j < n_trt; j ++) {
      idx = curr_vect[j];
      top += a[idx];
      lower_left += w[idx];
    }

    // compute the z-score
    out[k - start_pos] = top/sqrt(lower_left - lower_right);
  }

  return(out);
}

} // namespace inference
} // namespace sceptre
