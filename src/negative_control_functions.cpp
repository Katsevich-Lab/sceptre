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


// [[Rcpp::export]]
IntegerMatrix sample_combinations(int undercover_group_size, double n_pairs_to_sample, int N_NONZERO_TRT, int N_NONZERO_CNTRL, double n_possible_groups, IntegerMatrix n_nonzero_m, IntegerVector n_nonzero_tot) {
  // initialize set to store sampled vectors
  std::unordered_set<std::vector<int>, boost::hash<std::vector<int>>> H;
  std::vector<int> sample;

  // initialize random number generator
  std::mt19937 generator(4);
  std::uniform_real_distribution<double> distribution(0, 1);

  // initialze additional variables
  int B = 10000, n_nt_grnas = n_nonzero_m.nrow(), n_pairs_passed_qc = 0, n_nonzero_trt, n_nonzero_cntrl, gene_idx;
  double MULT_FACTOR = 5.0, n_genes = (double) n_nonzero_m.ncol();

  // 1. estimate the fraction of pairs that passes QC
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
List sample_undercover_pairs(IntegerMatrix n_nonzero_m, IntegerVector n_nonzero_tot, IntegerMatrix possible_groups_m, int n_pairs_to_sample, int N_NONZERO_TRT, int N_NONZERO_CNTRL) {
  // determine the number of grna groups, genes, and possible pairs
  int n_grna_groups = possible_groups_m.nrow();
  int undercover_grp_size = possible_groups_m.ncol();
  int n_genes = n_nonzero_m.ncol();
  int n_possible_pairs = n_grna_groups * n_genes;
  double n_possible_pairs_doub = (double) n_possible_pairs;
  int n_nonzero_trt, n_nonzero_cntrl;

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
    for (int j = 0; j < undercover_grp_size; j++) {
      n_nonzero_trt += n_nonzero_m(possible_groups_m(grna_group_idx, j), gene_idx);
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
