// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <algorithm>
#include <random>
#include <boost/container_hash/hash.hpp>
#include <unordered_set>
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


// output: a vector y whose length is equal to the number of cells, containing a true/false indicating presence/absence of a gene expression
void load_nonzero_posits(IntegerVector j, IntegerVector p, int column_idx, std::vector<bool>& y) {
  int start = p[column_idx];
  int end = p[column_idx + 1];
  for (int k = 0; k < y.size(); k ++) y[k] = false;
  for (int k = start; k < end; k ++) y[j[k]] = true;
  return;
}


// [[Rcpp::export]]
IntegerMatrix sample_combinations(int undercover_group_size, double n_pairs_to_sample, int N_NONZERO_TRT, int N_NONZERO_CNTRL, double n_possible_groups, IntegerMatrix n_nonzero_m, IntegerVector n_nonzero_tot, int N_POSSIBLE_GROUPS_THRESHOLD) {
  // initialize set to store sampled vectors
  std::unordered_set<std::vector<int>, boost::hash<std::vector<int>>> H;
  std::vector<int> sample;

  // initialize random number generator
  std::mt19937 generator(4);
  std::uniform_real_distribution<double> distribution(0, 1);

  // initialize additional variables
  int B = 10000, n_nt_grnas = n_nonzero_m.nrow(), n_pairs_passed_qc = 0, n_nonzero_trt, n_nonzero_cntrl, gene_idx;
  double MULT_FACTOR = 5.0, n_genes = (double) n_nonzero_m.ncol();

  // 1. estimate the fraction of negative control pairs that passes QC
  for (int i = 0; i < B; i ++) {
    sample = draw_wor_sample(n_nt_grnas, undercover_group_size, generator, distribution);
    gene_idx = floor(distribution(generator) * n_genes);
    n_nonzero_trt = 0;
    for (int j = 0; j < undercover_group_size; j++) {
      n_nonzero_trt += n_nonzero_m(sample[j], gene_idx);
    }
    n_nonzero_cntrl = n_nonzero_tot[gene_idx] - n_nonzero_trt;
    if ((n_nonzero_trt >= N_NONZERO_TRT) & (n_nonzero_cntrl >= N_NONZERO_CNTRL)) n_pairs_passed_qc ++;
  }
  double p_hat = n_pairs_passed_qc / ((double) B);

  // 2. estimate the number of grna groups to generate
  int n_grna_groups = ceil((MULT_FACTOR * n_pairs_to_sample)/(p_hat * n_genes));
  n_grna_groups = std::max(n_grna_groups, N_POSSIBLE_GROUPS_THRESHOLD);
  if ((double) n_grna_groups >= n_possible_groups) n_grna_groups = (int) n_possible_groups;

  // 3. sample n_grna_groups grna groups via rejection sampling
  while (H.size() < n_grna_groups) {
    sample = draw_wor_sample(n_nt_grnas, undercover_group_size, generator, distribution);
    std::sort(sample.begin(), sample.end(), std::less<int>());
    if (H.count(sample) == 0) H.insert(sample);
  }

  // 4. insert the sampled elements into output matrix
  IntegerMatrix out(n_grna_groups, undercover_group_size);
  int counter = 0;
  for (auto itr : H) {
    for (int j = 0; j < undercover_group_size; j ++) {
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
List sample_undercover_pairs(IntegerMatrix n_nonzero_m, IntegerVector n_nonzero_tot, IntegerMatrix possible_groups_m, int n_pairs_to_sample, int N_NONZERO_TRT, int N_NONZERO_CNTRL, bool low_moi, IntegerVector j, IntegerVector p, int n_cells, int n_genes, List indiv_nt_grna_idxs) {
  // determine the number of grna groups, genes, and possible pairs
  int n_grna_groups = possible_groups_m.nrow();
  int undercover_grp_size = possible_groups_m.ncol();
  int n_possible_pairs = n_grna_groups * n_genes;
  double n_possible_pairs_doub = (double) n_possible_pairs;
  int n_nonzero_trt, n_nonzero_cntrl;
  std::vector<bool> y(n_cells);
  IntegerVector curr_idxs;

  // initialize the vector of sampled elements
  std::vector<int> gene_idxs, grna_group_idxs, n_nonzero_trt_v, n_nonzero_cntrl_v;

  // initialize the vector of elements that could be sampled
  std::vector<int> x(n_possible_pairs);
  for (int i = 0; i < n_possible_pairs; i ++) x[i] = i;

  // initialize random number generator and variables related to random number generation
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

    if (low_moi) {
      // sum over elements of n nonzero mat
      for (int k = 0; k < undercover_grp_size; k++) {
        n_nonzero_trt += n_nonzero_m(possible_groups_m(grna_group_idx, k), gene_idx);
      }
    } else {
      // load vector of positions of nonzero entries into y
      load_nonzero_posits(j, p, gene_idx, y);
      // loop over the grnas in the group
      for (int k = 0; k < undercover_grp_size; k ++) {
        curr_idxs = indiv_nt_grna_idxs[possible_groups_m(grna_group_idx, k)];
        for (int l = 0; l < curr_idxs.size(); l++) {
          if (y[curr_idxs[l] - 1]) {
            n_nonzero_trt ++;
            y[curr_idxs[l] - 1] = false;
          }
        }
      }
    }
    n_nonzero_cntrl = n_nonzero_tot[gene_idx] - n_nonzero_trt;

    // if n_nonzero_trt and n_nonzero_cntrl exceed the thresholds, append the grna_group_idx and gene_idx
    if (n_nonzero_trt >= N_NONZERO_TRT && n_nonzero_cntrl >= N_NONZERO_CNTRL) {
      gene_idxs.push_back(gene_idx + 1);
      grna_group_idxs.push_back(grna_group_idx + 1);
      n_nonzero_trt_v.push_back(n_nonzero_trt);
      n_nonzero_cntrl_v.push_back(n_nonzero_cntrl);
    }

    // break out of loop if n_pairs_to_sample sampled
    if (gene_idxs.size() >= n_pairs_to_sample) break;
  }

  return(List::create(Named("response_idxs") = gene_idxs,
                      Named("grna_group_idxs") = grna_group_idxs,
                      Named("n_nonzero_trt_v") = n_nonzero_trt_v,
                      Named("n_nonzero_cntrl_v") = n_nonzero_cntrl_v));
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


// this function outputs four pieces:
// 1. N nonzero mat; this is the number of nonzero cells for each gene-NT gRNA pair (same regardless of control group)
// 2. n_nonzero_tot; this is the vector giving the number of nonzero cells per gene; when using the NT cells as the complement, we restrict our attention to the NT cells; when using the complement set as the control group, by contrast, we sum over all cells.
// 3,4. n nonzero trt, n nonzero cntrl vectors
// [[Rcpp::export]]
List compute_nt_nonzero_matrix_and_n_ok_pairs_v2(IntegerVector j, IntegerVector p, int n_cells, List grna_group_idxs, List indiv_nt_grna_idxs, IntegerVector all_nt_idxs, IntegerVector to_analyze_response_idxs, IntegerVector to_analyze_grna_idxs, bool compute_n_ok_pairs, bool control_group_complement) {
  // 0. initialize variables and objects
  std::vector<bool> y(n_cells);
  int n_nt_grnas = indiv_nt_grna_idxs.size(), n_genes = p.size() - 1, n_pairs = to_analyze_response_idxs.size(), pair_pointer = 0, n_nonzero = 0, curr_n_nonzero_cntrl = 0, curr_n_nonzero_trt = 0;
  IntegerVector curr_idxs;
  IntegerVector n_nonzero_tot(n_genes);
  bool n_cntrl_cells_ok;
  IntegerMatrix M(n_nt_grnas, n_genes);
  IntegerVector n_nonzero_trt(n_pairs), n_nonzero_cntrl(n_pairs);

  // 1. iterate over genes
  for (int column_idx = 0; column_idx < n_genes; column_idx++) {
    // load nonzero positions into the boolean vector y
    load_nonzero_posits(j, p, column_idx, y);
    // 1.2 iterate over nt grnas, adding n nonzero trt to M matrix
    for (int grna_idx = 0; grna_idx < n_nt_grnas; grna_idx ++) {
      n_nonzero = 0;
      curr_idxs = indiv_nt_grna_idxs[grna_idx];
      for (int k = 0; k < curr_idxs.size(); k ++) {
        if (!control_group_complement) { // NT cells control group
          if (y[all_nt_idxs[curr_idxs[k] - 1] - 1]) n_nonzero ++; // redirection
        } else { // discovery cells control group
          if (y[curr_idxs[k] - 1]) n_nonzero ++; // no redirection
        }
      }
      M(grna_idx, column_idx) = n_nonzero;
    }

    // 1.2 update n_nonzero_tot for this gene
    if (!control_group_complement) {
      n_nonzero_tot[column_idx] = 0;
      for (int k = 0; k < n_nt_grnas; k ++) n_nonzero_tot[column_idx] += M(k, column_idx);
    } else {
      n_nonzero_tot[column_idx] = p[column_idx + 1] - p[column_idx]; // or, the number of trues in the y vector
    }

    if (compute_n_ok_pairs) {
      // iterate through the discovery pairs containing this gene
      while (to_analyze_response_idxs[pair_pointer] - 1 == column_idx) {
        curr_idxs = grna_group_idxs[to_analyze_grna_idxs[pair_pointer] - 1];
        curr_n_nonzero_trt = 0;
        // first, get n nonzero trt
        for (int k = 0; k < curr_idxs.size(); k ++) if (y[curr_idxs[k] - 1]) curr_n_nonzero_trt ++;
        // second, get n nonzero cntrl
        if (!control_group_complement) { // NT cells
          curr_n_nonzero_cntrl = n_nonzero_tot[column_idx];
        } else { // complement set
          curr_n_nonzero_cntrl = n_nonzero_tot[column_idx] - curr_n_nonzero_trt;
        }
        n_nonzero_trt[pair_pointer] = curr_n_nonzero_trt;
        n_nonzero_cntrl[pair_pointer] = curr_n_nonzero_cntrl;
        pair_pointer ++;
      }
    }
  }

  return List::create(Named("n_nonzero_mat") = M, Named("n_nonzero_tot") = n_nonzero_tot,
                      Named("n_nonzero_trt") = n_nonzero_trt, Named("n_nonzero_cntrl") = n_nonzero_cntrl);
}
