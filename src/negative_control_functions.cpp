// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <algorithm>
#include <random>
#include <boost/container_hash/hash.hpp>
#include <unordered_set>
#include "shared_low_level_functions.h"
using namespace Rcpp;


std::vector<int> draw_wor_sample(int n, int k, std::mt19937& generator, std::uniform_real_distribution<double>& distribution) {
  std::vector<int> x(n);
  for (int i = 0; i < n; i ++) x[i] = i;
  std::vector<int> out(k);
  double n_doub = (double) n, u;
  int r;
  for (int i = 0; i < k; i ++) {
    u = distribution(generator);
    r = floor((n - (double) i) * u);
    out[i] = x[r];
    x[r] = x[n - i - 1];
  }
  return out;
}


// [[Rcpp::export]]
IntegerMatrix sample_combinations_v2(int calibration_group_size, double n_calibration_pairs, double n_possible_groups,
                                     int n_nt_grnas, double n_genes, int N_POSSIBLE_GROUPS_THRESHOLD, double p_hat) {
  // 0. initialize set to store sampled vectors
  std::unordered_set<std::vector<int>, boost::hash<std::vector<int>>> H;
  std::vector<int> sample;
  double MULT_FACTOR = 5.0;

  // 1. initialize random number generator
  std::mt19937 generator(4);
  std::uniform_real_distribution<double> distribution(0, 1);

  // 2. estimate the number of grna groups to generate
  int n_grna_groups = ceil((MULT_FACTOR * n_calibration_pairs)/(p_hat * n_genes));
  n_grna_groups = std::max(n_grna_groups, N_POSSIBLE_GROUPS_THRESHOLD);
  if ((double) n_grna_groups >= n_possible_groups) n_grna_groups = (int) n_possible_groups;

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


// [[Rcpp::export]]
IntegerMatrix iterate_over_combinations(int n_nt_grnas, int undercover_group_size, int n_possible_groups) {
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


// [[Rcpp::export]]
void increment_matrix(IntegerMatrix m) {
  int n_row = m.nrow();
  int n_col = m.ncol();
  for (int i = 0; i < n_row; i ++) {
    for (int j = 0; j < n_col; j ++) {
      m(i, j) ++;
    }
  }
  return;
}


// [[Rcpp::export]]
List sample_undercover_pairs_v2(IntegerMatrix n_nonzero_m, IntegerVector n_nonzero_tot, IntegerMatrix possible_groups_m, int n_genes,
                                int n_calibration_pairs, int n_nonzero_trt_thresh, int n_nonzero_cntrl_thresh, bool calculate_ess_using_m_matrix,
                                IntegerVector j, IntegerVector p, int n_cells_orig, int n_cells_sub, List indiv_nt_grna_idxs, IntegerVector cells_in_use) {
  // 0.1 define a few variables
  int n_grna_groups = possible_groups_m.nrow(), undercover_grp_size = possible_groups_m.ncol(), n_nonzero_trt, n_nonzero_cntrl;
  int n_possible_pairs = n_grna_groups * n_genes;
  double n_possible_pairs_doub = (double) n_possible_pairs;

  // 0.2 define variables related to the expression matrix
  std::vector<bool> y_sub(n_cells_sub), y_orig(n_cells_orig);
  std::vector<int> cells_in_use_zero_idx(n_cells_sub);
  IntegerVector curr_idxs;
  for (int i = 0; i < cells_in_use_zero_idx.size(); i ++) cells_in_use_zero_idx[i] = cells_in_use[i] - 1;

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
    n_nonzero_trt = 0;

    if (calculate_ess_using_m_matrix) {
      // sum over elements of n nonzero mat
      for (int k = 0; k < undercover_grp_size; k++) {
        n_nonzero_trt += n_nonzero_m(possible_groups_m(grna_group_idx, k), gene_idx);
      }
    } else {
      // load vector of positions of nonzero entries into y
      load_nonzero_posits(j, p, gene_idx, y_orig, y_sub, cells_in_use_zero_idx);
      // loop over the grnas in the group
      for (int k = 0; k < undercover_grp_size; k ++) {
        curr_idxs = indiv_nt_grna_idxs[possible_groups_m(grna_group_idx, k)];
        for (int l = 0; l < curr_idxs.size(); l++) {
          if (y_sub[curr_idxs[l] - 1]) {
            n_nonzero_trt ++;
            y_sub[curr_idxs[l] - 1] = false;
          }
        }
      }
    }
    n_nonzero_cntrl = n_nonzero_tot[gene_idx] - n_nonzero_trt;

    // if n_nonzero_trt and n_nonzero_cntrl exceed the thresholds, append the grna_group_idx and gene_idx
    if (n_nonzero_trt >= n_nonzero_trt_thresh && n_nonzero_cntrl >= n_nonzero_cntrl_thresh) {
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
