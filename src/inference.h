#ifndef SCEPTRE_INFERENCE_H
#define SCEPTRE_INFERENCE_H

#include <Rcpp.h>
#include <array>
#include <vector>

// Numerical pieces of the full-score association test. Inside C++, observed and
// resampled cell positions start at 0. The R-facing run_association_test()
// converts observed R positions and reads the draws through their external pointer.
namespace sceptre {
namespace inference {

// A shared collection of B random draws, each listing selected cell positions.
// Permutation tests use the first n_trt positions; CRT uses the whole list,
// whose length can vary. A hybrid draw can include later positions intended for
// a larger fit, so a smaller fit must read only its own prefix. Reuse draws
// unchanged across genes and stages. Entries are not original matrix cell IDs.
typedef std::vector<std::vector<int>> DrawBank;

// Check that score arrays describe the same cells. Callers check observed cell
// positions once. Samplers make the random positions valid; permutation tests
// read only the requested prefix, and CRT reads the complete selected draw.
void check_score_dimensions(const Rcpp::NumericVector& a,
                            const Rcpp::NumericVector& w,
                            const Rcpp::NumericMatrix& D);

// Copy observed R positions and subtract one, keeping their order and repeats.
// Reject missing/nonpositive values; the R entry checks positions against the
// current fit's cell count before calling the numerical test.
Rcpp::IntegerVector observed_positions_zero_based(const Rcpp::IntegerVector& positions);

// Fold-change estimate and residual-based standard error, both before taking logs.
struct EffectEstimate {
  double fold_change;
  double standard_error;
};

struct TailResult {
  double p;
  int stage;
  std::array<double, 3> sn_params;
  // Null scores from the final stage used, in the order left by tail fitting.
  // This is a separate vector; sorting it does not change the shared draws.
  std::vector<double> null_statistics;
};

// Estimate effect from counts y and means mu in the same cell order, using n_trt
// treated positions starting at 0. Leave inputs unchanged and preserve sum order
// and ordinary NaN/Inf results for undefined calculations.
EffectEstimate estimate_effect(const Rcpp::NumericVector& y,
                               const Rcpp::NumericVector& mu,
                               const Rcpp::IntegerVector& trt_idxs, int n_trt);

// Calculate the standardized score from one a/w value per cell, D with one column
// per cell, and treated positions starting at 0. Inputs are already prepared:
// this function does not fit models, generate draws or repair degenerate values.
double observed_score(const Rcpp::NumericVector& a,
                      const Rcpp::NumericVector& w,
                      const Rcpp::NumericMatrix& D,
                      const Rcpp::IntegerVector& trt_idxs);

// Calculate B null scores from draws start_pos through start_pos+B-1.
// Use n_trt positions per draw unless use_all_cells requests each complete CRT
// draw. The caller checks score-array dimensions before allocating C++ results.
std::vector<double> null_score_segment(const Rcpp::NumericVector& a,
                                       const Rcpp::NumericVector& w,
                                       const Rcpp::NumericMatrix& D,
                                       int start_pos, int B, int n_trt,
                                       bool use_all_cells,
                                       const DrawBank& synthetic_idxs);

// Count ties in empirical tails and add one to tail counts and total draw count.
// For a two-sided test, double the smaller tail and cap at 1.
double compute_empirical_p_value(const std::vector<double>& null_statistics,
                                 double z_orig, int side);

// Fit a skew-normal tail to a separate vector of at least 100 scores. Sort that
// vector, possibly leaving it descending after a left-tail check. Return
// xi/omega/alpha/p, with p = -1 requesting empirical fallback. Full draws stay
// unchanged.
std::vector<double> fit_and_evaluate_skew_normal(double z_orig,
                                               std::vector<double>& null_statistics,
                                               int side_code);

// Use the staged test with cutoff .02, an optional skew-normal fit and empirical
// fallback. When B3 = 0, keep the previous stage's scores rather than joining
// scores from multiple stages.
TailResult evaluate_staged_tail(double z_orig,
                                const Rcpp::NumericVector& a,
                                const Rcpp::NumericVector& w,
                                const Rcpp::NumericMatrix& D,
                                int n_trt, bool use_all_cells,
                                const DrawBank& synthetic_idxs,
                                int B1, int B2, int B3,
                                bool fit_parametric_curve, int side_code);

// Copy the finished C++ tail result into R and free it, including if R errors.
// Before computing it, the caller must have PROTECTed an R continuation token.
// After releasing its C++ result pointer, call this directly: any intervening R
// allocation/error could skip C++ cleanup and leak the result. Conversion uses
// the raw R API with registered cleanup. UNPROTECT the token once after return;
// R restores its protection stack itself on the protected error path.
SEXP association_result_record(const EffectEstimate& effect, double z_orig,
                               TailResult* tail, bool return_resampling_dist,
                               SEXP continuation);

} // namespace inference
} // namespace sceptre

#endif
