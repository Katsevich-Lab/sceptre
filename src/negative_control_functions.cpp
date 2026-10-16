// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <algorithm>
#include <random>
#include <boost/container_hash/hash.hpp>
#include <unordered_set>
#include "pair_eligibility.h"
#include "sampling_safety.h"
using namespace Rcpp;

// Calibration uses private std engines and consumes no R random draws. The
// generated Rcpp RNGScope is retained, including initialization of an absent
// R seed; an existing R stream is not advanced by these native draws.

namespace {
// Calibration needs nonempty NT guides and 1<=group_size<=active NT count.
// Guard enumeration/sampling before indexing or allocating; no metadata fallback.
void check_group_size(int n_nt_grnas, int group_size) {
  if (n_nt_grnas < 1 || group_size < 1 || group_size > n_nt_grnas)
    stop("Calibration group size must be between 1 and the active NT group count.");
}

// Double-valued count guard before narrowing to native int; finite, integral,
// positive and within capacity. This does not choose a smaller feasible request.
void positive_integer(double value, const char* name) {
  if (!std::isfinite(value) || value < 1 || value != std::floor(value) ||
      value > std::numeric_limits<int>::max())
    stop("%s must be a positive integer within native capacity.", name);
}
}


// Distinct zero-based NT-guide subset in draw order, using caller-owned std
// mt19937/uniform state and O(n) scratch. Caller establishes 0<=k<=n. This std
// stream belongs to calibration; it is deliberately not the Boost bank stream.
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


// Sample a bounded set of distinct unordered NT-guide combinations, returning
// group x member integer matrix with zero-based active-NT-guide positions.
// Group count is min(choose(n,k), max(threshold,ceil(5*pairs/(p_hat*genes))));
// validate empirical p_hat and capacity before narrowing/allocation. Each subset
// is sorted for uniqueness; result rows follow the existing unordered_set
// iteration order, not a newly imposed order. Fixed std mt19937(4), no R draws
// or input mutation. Preserve rejection/draw order and the existing heuristic.
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
  // 0. initialize set to store sampled vectors
  std::unordered_set<std::vector<int>, boost::hash<std::vector<int>>> H;
  std::vector<int> sample;
  double MULT_FACTOR = 5.0;

  // 1. initialize random number generator
  std::mt19937 generator(4);
  std::uniform_real_distribution<double> distribution(0, 1);

  // 2. estimate the number of grna groups to generate
  double selected_groups = ceil((MULT_FACTOR * n_calibration_pairs)/(p_hat * n_genes));
  selected_groups = std::max(selected_groups, static_cast<double>(N_POSSIBLE_GROUPS_THRESHOLD));
  selected_groups = std::min(selected_groups, n_possible_groups);
  positive_integer(selected_groups, "Selected calibration group count");
  const int n_grna_groups = static_cast<int>(selected_groups);
  sceptre::sampling::matrix_extent(n_grna_groups, calibration_group_size);

  // 3. sample n_grna_groups grna groups via rejection sampling
  while (H.size() < n_grna_groups) {
    sample = draw_wor_sample(n_nt_grnas, calibration_group_size, generator, distribution);
    std::sort(sample.begin(), sample.end(), std::less<int>());
    if (H.count(sample) == 0) H.insert(sample);
  }

  // 4. insert the sampled elements into output matrix
  IntegerMatrix out(n_grna_groups, calibration_group_size);
  int counter = 0;
  for (auto itr : H) {
    for (int j = 0; j < calibration_group_size; j ++) {
      out(counter, j) = itr[j];
    }
    counter ++;
  }
  return out;
}


// Enumerate all size-k subsets of active NT guides into a group x member
// zero-based integer matrix. Preserve the current next_permutation visitation
// order (including its first step/wrap), rather than sorting output rows.
// Check choose(n,k), dimensions and capacity before filling. No RNG/mutation.
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


// Convert an explicitly owned scratch combination matrix from zero-based NT
// guide positions to R one-based label positions. MUTATES the supplied matrix,
// preserving shape/member/row order; the R calibration constructor calls this
// only after sampling has finished with the zero-based representation.
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


// Sample QC-eligible group/response pairs without replacement. Count matrix
// is active-NT-group x all-response; possible_groups_m is candidate-group x
// member, with ZERO-based row positions in that matrix. Summing member counts
// is correct for low-MOI/disjoint groups or the caller's precomputed union row;
// high-MOI overlap must be unioned before entry. Totals share response columns.
// Enumerate candidates by group*n_genes+response, traverse a fixed std seed-4
// partial Fisher-Yates stream and retain shared legacy_pair_eligible outcomes.
// Return one-based response/candidate-group positions plus counts in acceptance
// order; may exhaust before requested pairs. Widen/check products and sums
// before allocation/narrowing. Inputs read only, thresholds inclusive, no R RNG;
// do not replace this stream/eligibility policy with a new acceptance heuristic.
// [[Rcpp::export]]
List sample_calibration_pairs(IntegerMatrix n_nonzero_m, IntegerVector n_nonzero_tot, IntegerMatrix possible_groups_m,
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

  // 0.3 initialize the vector of sampled elements and elements that could be sampld
  std::vector<int> gene_idxs, grna_group_idxs, n_nonzero_trt_v, n_nonzero_cntrl_v;
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
      treatment_count += n_nonzero_m(possible_groups_m(grna_group_idx, k), gene_idx);
      if (treatment_count > std::numeric_limits<int>::max())
        stop("Calibration treatment count exceeds native integer capacity.");
    }
    if (treatment_count > n_nonzero_tot[gene_idx])
      stop("Calibration treatment count exceeds its total nonzero count.");
    n_nonzero_trt = static_cast<int>(treatment_count);
    n_nonzero_cntrl = n_nonzero_tot[gene_idx] - n_nonzero_trt;

    // if n_nonzero_trt and n_nonzero_cntrl exceed the thresholds, append the grna_group_idx and gene_idx
    if (legacy_pair_eligible(n_nonzero_trt, n_nonzero_cntrl, n_nonzero_trt_thresh, n_nonzero_cntrl_thresh)) {
      gene_idxs.push_back(gene_idx + 1);
      grna_group_idxs.push_back(grna_group_idx + 1);
      n_nonzero_trt_v.push_back(n_nonzero_trt);
      n_nonzero_cntrl_v.push_back(n_nonzero_cntrl);
    }

    // break out of loop if n_calibration_pairs sampled
    if (gene_idxs.size() >= n_calibration_pairs) break;
  }

  return(List::create(Named("response_idxs") = gene_idxs,
                      Named("grna_group_idxs") = grna_group_idxs,
                      Named("n_nonzero_trt_v") = n_nonzero_trt_v,
                      Named("n_nonzero_cntrl_v") = n_nonzero_cntrl_v));
}
