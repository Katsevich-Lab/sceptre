// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <algorithm>
#include <random>
#include <boost/container_hash/hash.hpp>
#include <unordered_set>
#include "pair_eligibility.h"
#include "sampling_safety.h"
#include "result_lifetime.h"
using namespace Rcpp;

// Calibration uses its own standard-library random generators, separate from R.
// The generated Rcpp wrapper can initialize an absent R seed, but these C++
// draws do not advance an R random stream that already exists.

namespace {
// Calibration needs NT guides with remaining cells, and a group size between 1
// and that active guide count. Check before indexing or allocating combinations;
// metadata rows alone do not show which guides still have cells after QC.
void check_group_size(int n_nt_grnas, int group_size) {
  if (n_nt_grnas < 1 || group_size < 1 || group_size > n_nt_grnas)
    stop("Calibration group size must be between 1 and the active NT group count.");
}

// Check a double-valued request before converting to int: it must be finite,
// positive, whole and within capacity. Error rather than choose a smaller request.
void positive_integer(double value, const char* name) {
  if (!std::isfinite(value) || value < 1 || value != std::floor(value) ||
      value > std::numeric_limits<int>::max())
    stop("%s must be a positive integer within native capacity.", name);
}

// Store only pairs that actually pass QC. Reserving the entire request size can
// waste substantial memory when most candidates fail. The C++ result holds these
// vectors until protected conversion allocates R vectors of the accepted size.
struct CalibrationPairResult {
  std::vector<int> response_idxs, grna_group_idxs, n_nonzero_trt_v, n_nonzero_cntrl_v;
};

// Copy accepted pairs into R using the raw R API. The caller has already set up
// R_UnwindProtect with a protected continuation token, so an R allocation error
// still frees the C++ result. Do not create Rcpp cleanup objects or C++ buffers
// here: an R error could skip their destructors.
SEXP calibration_pair_record(void* data) {
  const CalibrationPairResult* result = static_cast<CalibrationPairResult*>(data);
  const R_xlen_t size = result->response_idxs.size();
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
  const char* labels[] = {"response_idxs", "grna_group_idxs", "n_nonzero_trt_v", "n_nonzero_cntrl_v"};
  const std::vector<int>* values[] = {&result->response_idxs, &result->grna_group_idxs,
                                    &result->n_nonzero_trt_v, &result->n_nonzero_cntrl_v};
  for (int field = 0; field < 4; ++field) {
    SET_STRING_ELT(names, field, Rf_mkChar(labels[field]));
    SET_VECTOR_ELT(out, field, Rf_allocVector(INTSXP, size));
    int* target = INTEGER(VECTOR_ELT(out, field));
    for (R_xlen_t i = 0; i < size; ++i) target[i] = (*values[field])[i];
  }
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

// Free the C++ result once, after conversion or before R continues an error.
// Deleting vectors must not allocate or throw. Between releasing the C++ owner
// and entering R_UnwindProtect, do not make an R call that could error and leak it.
void delete_calibration_pair_result(void* data, Rboolean) noexcept {
  delete static_cast<CalibrationPairResult*>(data);
}
}


// Sample k distinct NT-guide positions among n guides, starting at 0, in draw
// order. Reuse the caller's std::mt19937 generator and O(n) working array.
// The caller checks 0 <= k <= n. Calibration deliberately uses this standard-
// library generator, rather than the Boost generator used for association draws.
std::vector<int> draw_wor_sample(int n, int k, std::mt19937& generator, std::uniform_real_distribution<double>& distribution) {
  std::vector<int> x(n);
  for (int i = 0; i < n; i ++) x[i] = i;
  std::vector<int> out(k);
  double u;
  int r;
  for (int i = 0; i < k; i ++) {
    u = distribution(generator);
    r = floor((n - (double) i) * u);
    out[i] = x[r];
    x[r] = x[n - i - 1];
  }
  return out;
}


// Sample distinct combinations of k active NT guides. Return a matrix with groups
// in rows and members in columns, as guide positions starting at 0. Choose up to
// min(choose(n,k), max(threshold, ceil(5*pairs/(p_hat*genes)))) groups, following
// the current heuristic. Check p_hat and sizes before conversion/allocation.
// Sort each subset to detect repeated combinations, but keep the existing
// unordered_set iteration order for output rows. Start std::mt19937 at seed 4
// and preserve rejection/draw order; R's random stream and inputs are unchanged.
// [[Rcpp::export]]
IntegerMatrix sample_calibration_groups(int calibration_group_size, double n_calibration_pairs, double n_possible_groups,
                                     int n_nt_grnas, double n_genes, int N_POSSIBLE_GROUPS_THRESHOLD, double p_hat) {
  check_group_size(n_nt_grnas, calibration_group_size);
  positive_integer(n_calibration_pairs, "n_calibration_pairs");
  positive_integer(n_genes, "n_genes");
  positive_integer(N_POSSIBLE_GROUPS_THRESHOLD, "Group threshold");
  if (!std::isfinite(p_hat) || p_hat <= 0 || p_hat > 1)
    stop("Cannot sample calibration groups: the empirical candidate rate p_hat must be finite and in (0, 1].");
  const double represented_groups = R::choose(n_nt_grnas, calibration_group_size);
  if (std::isnan(n_possible_groups) || n_possible_groups < 1 ||
      n_possible_groups != represented_groups)
    stop("n_possible_groups must equal the represented combination count.");
  double MULT_FACTOR = 5.0;

  // Check the output size before creating C++ containers that need destruction.
  double selected_groups = ceil((MULT_FACTOR * n_calibration_pairs)/(p_hat * n_genes));
  selected_groups = std::max(selected_groups, static_cast<double>(N_POSSIBLE_GROUPS_THRESHOLD));
  selected_groups = std::min(selected_groups, n_possible_groups);
  positive_integer(selected_groups, "Selected calibration group count");
  const int n_grna_groups = static_cast<int>(selected_groups);
  sceptre::sampling::matrix_extent(n_grna_groups, calibration_group_size);

  IntegerMatrix out(n_grna_groups, calibration_group_size);
  int* output = INTEGER(out);
  {
    // All R allocations and data loading are finished. Use only C++ while the set
    // and working vectors are alive, so a C++ exception frees them normally and an
    // R error cannot skip their destructors.
    std::unordered_set<std::vector<int>, boost::hash<std::vector<int>>> H;
    std::vector<int> sample;
    std::mt19937 generator(4);
    std::uniform_real_distribution<double> distribution(0, 1);

    while (H.size() < n_grna_groups) {
      sample = draw_wor_sample(n_nt_grnas, calibration_group_size, generator, distribution);
      std::sort(sample.begin(), sample.end(), std::less<int>());
      if (H.count(sample) == 0) H.insert(sample);
    }

    int counter = 0;
    for (auto itr : H) {
      for (int j = 0; j < calibration_group_size; j ++) {
        output[counter + static_cast<R_xlen_t>(n_grna_groups) * j] = itr[j];
      }
      counter ++;
    }
  }
  return out;
}


// List all combinations of k active NT guides, with groups in rows and members
// in columns, as positions starting at 0. Keep the current next_permutation
// order, including its first step and wraparound, rather than sorting the rows.
// Check choose(n,k) and dimensions before filling the output. No random draws.
// [[Rcpp::export]]
IntegerMatrix iterate_over_combinations(int n_nt_grnas, int undercover_group_size, int n_possible_groups) {
  check_group_size(n_nt_grnas, undercover_group_size);
  sceptre::sampling::nonnegative(n_possible_groups, "n_possible_groups");
  if (static_cast<double>(n_possible_groups) != R::choose(n_nt_grnas, undercover_group_size))
    stop("n_possible_groups must equal the represented combination count.");
  sceptre::sampling::matrix_extent(n_possible_groups, undercover_group_size);
  // initialize the IntegerMatrix out
  IntegerMatrix m(n_possible_groups, undercover_group_size);

  // initialize the boolean indicator vector
  std::vector<int> x(n_nt_grnas, false);
  for (int i = 0; i < undercover_group_size; i ++) x[i] = true;
  int counter;

  // fill the matrix m with subsets of size undercover_group_size
  for (int i = 0; i < n_possible_groups; i ++) {
    std::next_permutation(x.begin(), x.end());
    counter = 0;
    for (int j = 0; j < n_nt_grnas; j++) if (x[j]) m(i, counter++) = j;
  }
  return(m);
}


// Add one to the supplied working combination matrix so its entries can index
// R guide labels. This changes the input matrix in place, preserving row/member
// order and shape. The R caller must supply its own temporary matrix and call
// this only after finishing all sampling that uses positions starting at 0.
// [[Rcpp::export]]
void increment_group_indices_in_place(IntegerMatrix m) {
  int n_row = m.nrow();
  int n_col = m.ncol();
  for (int i = 0; i < n_row; i ++) {
    for (int j = 0; j < n_col; j ++) {
      m(i, j) ++;
    }
  }
  return;
}


// Sample distinct calibration guide-group/gene pairs that pass the shared QC rule.
// Count-matrix rows are active NT groups and columns are genes. Each row of
// possible_groups_m lists positions in those count rows, starting at 0.
// Summing member counts is valid for disjoint low-MOI groups, or for a precomputed
// union row. For overlapping high-MOI guides, the caller must first count their
// union so cells are not counted twice. Totals use the same gene-column order.
// Encode candidates as group*n_genes + gene and visit them with the current
// partial Fisher-Yates stream, using std seed 4. Keep the shared inclusive QC
// thresholds and this sampling order. Return accepted gene/group positions
// starting at 1 and their counts, in acceptance order. Candidates may run out
// before the requested size. Check products and sums before allocating or
// converting to int. Caller inputs and R's random stream stay unchanged.
// [[Rcpp::export]]
SEXP sample_calibration_pairs(IntegerMatrix n_nonzero_m, IntegerVector n_nonzero_tot, IntegerMatrix possible_groups_m,
                                int n_genes, int n_calibration_pairs, int n_nonzero_trt_thresh, int n_nonzero_cntrl_thresh) {
  positive_integer(n_genes, "n_genes");
  positive_integer(n_calibration_pairs, "n_calibration_pairs");
  sceptre::sampling::nonnegative(n_nonzero_trt_thresh, "Treatment count threshold");
  sceptre::sampling::nonnegative(n_nonzero_cntrl_thresh, "Control count threshold");
  if (n_nonzero_m.ncol() != n_genes || n_nonzero_tot.size() != n_genes)
    stop("Calibration gene counts and totals must have matching gene dimensions.");
  if (possible_groups_m.nrow() < 1 || possible_groups_m.ncol() < 1)
    stop("Calibration group definitions must be nonempty.");
  for (int row : possible_groups_m) {
    if (row < 0 || row >= n_nonzero_m.nrow())
      stop("Calibration group row is outside the zero-based count matrix range.");
  }
  for (int count : n_nonzero_m) sceptre::sampling::nonnegative(count, "Group nonzero count");
  for (int count : n_nonzero_tot) sceptre::sampling::nonnegative(count, "Total nonzero count");
  // 0.1 define a few variables
  int n_grna_groups = possible_groups_m.nrow(), undercover_grp_size = possible_groups_m.ncol(), n_nonzero_trt, n_nonzero_cntrl;
  const int64_t pair_count = static_cast<int64_t>(n_grna_groups) * n_genes;
  if (pair_count > std::numeric_limits<int>::max() ||
      static_cast<uint64_t>(pair_count) > std::vector<int>().max_size())
    stop("Candidate pair count exceeds native integer allocation capacity.");
  int n_possible_pairs = static_cast<int>(pair_count);
  double n_possible_pairs_doub = (double) n_possible_pairs;

  // Load all R input data and protect the error-continuation token before creating
  // the C++ result; an R error must not skip cleanup of already allocated vectors.
  const int* counts = INTEGER(n_nonzero_m);
  const int* totals = INTEGER(n_nonzero_tot);
  const int* group_members = INTEGER(possible_groups_m);
  const int count_rows = n_nonzero_m.nrow();
  SEXP continuation;
  CalibrationPairResult* owned_result;
  {
    sceptre::ResultContinuation token;
    std::unique_ptr<CalibrationPairResult> result(new CalibrationPairResult);
    std::vector<int> x(n_possible_pairs);
    for (int i = 0; i < n_possible_pairs; i ++) x[i] = i;

    // 0.4 initialize random number generator and variables related to random number generation
    std::mt19937 generator(4);
    std::uniform_real_distribution<double> distribution(0, 1);
    double u;
    int r, idx, grna_group_idx, gene_idx;

    // sample pairs
    for (int i = 0; i < n_possible_pairs; i ++) {
      // sample a random number in the range [0, ..., n_possible_pairs - i - 1]
      u = distribution(generator);
      r = floor((n_possible_pairs_doub - (double) i) * u);
      idx = x[r];
      x[r] = x[n_possible_pairs - i - 1];

      // obtain the grna group index and gene index
      grna_group_idx = idx/n_genes;
      gene_idx = idx % n_genes;

      // check if the grna group-gene pair passes pairwise QC
      int64_t treatment_count = 0;
      for (int k = 0; k < undercover_grp_size; k++) {
        treatment_count += counts[group_members[grna_group_idx + static_cast<R_xlen_t>(n_grna_groups) * k] + static_cast<R_xlen_t>(count_rows) * gene_idx];
        if (treatment_count > std::numeric_limits<int>::max())
          stop("Calibration treatment count exceeds native integer capacity.");
      }
      if (treatment_count > totals[gene_idx])
        stop("Calibration treatment count exceeds its total nonzero count.");
      n_nonzero_trt = static_cast<int>(treatment_count);
      n_nonzero_cntrl = totals[gene_idx] - n_nonzero_trt;

      // Keep the pair if both counts meet their inclusive QC thresholds.
      if (legacy_pair_eligible(n_nonzero_trt, n_nonzero_cntrl, n_nonzero_trt_thresh, n_nonzero_cntrl_thresh)) {
        result->response_idxs.push_back(gene_idx + 1);
        result->grna_group_idxs.push_back(grna_group_idx + 1);
        result->n_nonzero_trt_v.push_back(n_nonzero_trt);
        result->n_nonzero_cntrl_v.push_back(n_nonzero_cntrl);
      }

      // break out of loop if n_calibration_pairs sampled
      if (result->response_idxs.size() >= n_calibration_pairs) break;
    }
    owned_result = result.release();
    continuation = token.release();
  }
  // Pass the released result directly to protected R conversion. No R call may
  // intervene before conversion takes responsibility for freeing its C++ memory.
  SEXP out = R_UnwindProtect(calibration_pair_record, owned_result,
                             delete_calibration_pair_result, owned_result, continuation);
  UNPROTECT(1);
  return out;
}
