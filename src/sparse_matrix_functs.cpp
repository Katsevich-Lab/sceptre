#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

void load_sparse_vector(IntegerVector l, IntegerVector p, NumericVector x, int idx, int dim, NumericVector out) {
  int start = p[idx], end = p[idx + 1];
  for (int k = 0; k < dim; k ++) out[k] = 0.0;
  for (int k = start; k < end; k ++) out[l[k]] = x[k];
}


// [[Rcpp::export]]
NumericVector load_csr_row(IntegerVector j, IntegerVector p, NumericVector x, int row_idx, int n_cells) {
  NumericVector out(n_cells);
  load_sparse_vector(j, p, x, row_idx - 1, n_cells, out);
  return out;
}


// [[Rcpp::export]]
IntegerVector obtain_pointer_vector(IntegerVector i, int dim) {
  IntegerVector p(dim + 1);
  p[0] = 0;
  // special case of dim = 1
  if (dim == 1) {
    p[1] = i.size();
  } else {
    int curr_idx = 0, counter = 0, j = 0;
    for (curr_idx = 0; curr_idx < dim; curr_idx ++) {
      while (i[j] == curr_idx) {
        counter ++;
        j ++;
      }
      p[curr_idx + 1] = counter;
    }
  }
  return p;
}


// [[Rcpp::export]]
List compute_cell_covariates_cpp(IntegerVector i, IntegerVector p, NumericVector x, int n_genes,
                                 int n_cells, IntegerVector mt_gene_idxs, bool compute_p_mito, bool compute_max_feature) {
  IntegerVector n_nonzero(n_cells), max_feature(n_cells);
  NumericVector n_umi(n_cells), p_mito(n_cells), y(n_genes), frac_umis_max_feature(n_cells);
  List out;

  int n_nonzero_count, curr_max_feature;
  double n_umi_count, n_umi_count_mito, umi_count_curr_max_feature;
  for (int k = 0; k < mt_gene_idxs.size(); k ++) mt_gene_idxs[k] --;

  for (int k = 0; k < n_cells; k ++) {
    curr_max_feature = 0;
    umi_count_curr_max_feature = 0;
    n_umi_count = 0;
    n_nonzero_count = 0;
    load_sparse_vector(i, p, x, k, n_genes, y);
    for (int r = 0; r < n_genes; r ++) {
      n_umi_count += y[r];
      if (y[r] > 0.5) n_nonzero_count ++;
    }

    // compute p mito
    if (compute_p_mito) {
      n_umi_count_mito = 0;
      for (int r = 0; r < mt_gene_idxs.size(); r ++) n_umi_count_mito += y[mt_gene_idxs[r]];
      p_mito[k] = n_umi_count_mito/n_umi_count;
    }

    // compute max feature, fraction expressed
    if (compute_max_feature) {
      curr_max_feature = 0;
      umi_count_curr_max_feature = 0;
      for (int r = 0; r < n_genes; r ++) {
        if (y[r] > umi_count_curr_max_feature) {
          curr_max_feature = r;
          umi_count_curr_max_feature = y[r];
        }
      }
      max_feature[k] = curr_max_feature;
      frac_umis_max_feature[k] = umi_count_curr_max_feature/n_umi_count;
    }

    n_nonzero[k] = n_nonzero_count;
    n_umi[k] = n_umi_count;
  }

  out = List::create(Named("n_nonzero") = n_nonzero, Named("n_umi") = n_umi, Named("p_mito") = p_mito,
                     Named("max_feature") = max_feature, Named("frac_umis_max_feature") = frac_umis_max_feature);
  return out;
}


// [[Rcpp::export]]
List compute_colwise_max(IntegerVector i, IntegerVector p, NumericVector x, int n_cells, NumericVector grna_lib_size) {
  IntegerVector assignment_vect(n_cells);
  NumericVector frac_umis(n_cells);
  int p_start, p_end, curr_maximizer;
  double curr_max;

  for (int k = 0; k < n_cells; k ++) {
    p_start = p[k];
    p_end = p[k+1];
    curr_max = -1.0;
    curr_maximizer = -1;
    for (int l = p_start; l < p_end; l ++) {
      if (x[l] > curr_max) {
        curr_max = x[l];
        curr_maximizer = i[l];
      }
    }
    assignment_vect[k] = curr_maximizer + 1;
    frac_umis[k] = curr_max/grna_lib_size[k];
  }
  return(List::create(Named("assignment_vect") = assignment_vect, Named("frac_umis") = frac_umis));
}


// [[Rcpp::export]]
IntegerVector compute_n_grnas_per_cell_vector(List grna_assignments, int n_cells) {
  IntegerVector out(n_cells), current_assignments;
  for (int i = 0; i < grna_assignments.size(); i ++) {
    current_assignments = grna_assignments[i];
    for (int j = 0; j < current_assignments.size(); j ++) {
      out[current_assignments[j] - 1] ++;
    }
  }
  return out;
}


// [[Rcpp::export]]
void increment_vector(IntegerVector x, int value) {
  for (int i = 0; i < x.size(); i++) x[i] += value;
  return;
}


// [[Rcpp::export]]
IntegerVector threshold_count_matrix(IntegerVector j, IntegerVector p, NumericVector x, int row_idx, double threshold) {
  IntegerVector out;
  int n_gte_threshold = 0, counter = 0, start = p[row_idx - 1], end = p[row_idx];
  for (int k = start; k < end; k ++) if (x[k] >= threshold) n_gte_threshold ++;
  out = IntegerVector(n_gte_threshold);
  for (int k = start; k < end; k ++) if (x[k] >= threshold) out[counter ++] = j[k] + 1;
  return out;
}
