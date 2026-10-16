#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// Load one zero-based compressed segment into caller-owned dense scratch out
// (length dim), resetting it first. Slots/pointers are canonical zero-based
// Matrix CSR/CSC coordinates; adapters establish valid layout and dimensions.
// Explicit zero/fractional VALUES are loaded; no threshold or stored-position
// counting occurs. Canonical segments have unique indices (otherwise last wins).
// Mutates only out, never source slots, and allocates nothing/no RNG.
void load_sparse_vector(IntegerVector l, IntegerVector p, NumericVector x, int idx, int dim, NumericVector out) {
  int start = p[idx], end = p[idx + 1];
  for (int k = 0; k < dim; k ++) out[k] = 0.0;
  for (int k = start; k < end; k ++) out[l[k]] = x[k];
}


// Dense numeric row in original source-cell order from valid CSR j/p/x.
// row_idx is R one-based, converted once for load_sparse_vector; n_cells is
// full source width. Returns a new vector, retaining loaded values and implicit
// zeros; storage selects/reorders requested cells afterward. No input mutation.
// [[Rcpp::export]]
NumericVector load_csr_row(IntegerVector j, IntegerVector p, NumericVector x, int row_idx, int n_cells) {
  NumericVector out(n_cells);
  load_sparse_vector(j, p, x, row_idx - 1, n_cells, out);
  return out;
}


// Build dim+1 zero-based compressed pointers from an already nondecreasing
// stream of zero-based major indices in [0,dim). Empty rows/columns get repeated
// offsets; dim=1 uses all entries. Import/storage callers establish sorted valid
// coordinates before entry. Repetitions remain stored entries, not coalesced.
// Returns a new vector, with no mutation, value inspection, RNG or file effect.
// [[Rcpp::export]]
IntegerVector obtain_pointer_vector(IntegerVector i, int dim) {
  int i_size = i.size();
  IntegerVector p(dim + 1);
  p[0] = 0;
  // special case of dim = 1
  if (dim == 1) {
    p[1] = i_size;
  } else {
    int curr_idx = 0, counter = 0;
    for (curr_idx = 0; curr_idx < dim; curr_idx ++) {
      while (counter < i_size && i[counter] == curr_idx) counter ++;
      p[curr_idx + 1] = counter;
    }
  }
  return p;
}


// Feature-by-source-cell CSC summaries: i/p/x are canonical zero-based slots,
// n_genes rows/n_cells columns; output vectors follow original cell order.
// Column summaries load values (not stored-position counts): nonzero means
// value > 0.5. UMIs sum all loaded values, including fractional counts; optional
// mito fractions use R one-based feature rows. Optional max-feature returns a
// zero-based first strictly maximal row; all-zero columns keep row 0/NaN ratio.
// Unrequested optional fields remain their allocated zeros. Active mitochondrial indices are copied from R one-based rows
// only after all R allocation/materialization, including the named result.
// The native compute scope uses cached pointers and ends before return. Source
// i/p/x and all R inputs remain unchanged; preserve ties and zero-total ratios.
// [[Rcpp::export]]
List compute_cell_covariates_cpp(IntegerVector i, IntegerVector p, NumericVector x, int n_genes,
                                 int n_cells, IntegerVector mt_gene_idxs, bool compute_p_mito, bool compute_max_feature) {
  IntegerVector n_nonzero(n_cells), max_feature(n_cells);
  NumericVector n_umi(n_cells), p_mito(n_cells), y(n_genes), frac_umis_max_feature(n_cells);
  List out = List::create(Named("n_nonzero") = n_nonzero, Named("n_umi") = n_umi, Named("p_mito") = p_mito,
                          Named("max_feature") = max_feature, Named("frac_umis_max_feature") = frac_umis_max_feature);
  const int* mt_rows = nullptr;
  R_xlen_t n_mt = 0;
  if (compute_p_mito) {
    mt_rows = INTEGER(mt_gene_idxs);
    n_mt = mt_gene_idxs.size();
    for (R_xlen_t row = 0; row < n_mt; ++row)
      if (mt_rows[row] < 1 || mt_rows[row] > n_genes)
        stop("Mitochondrial row is outside the one-based feature range.");
  }
  const int* indices = INTEGER(i);
  const int* pointers = INTEGER(p);
  const double* values = REAL(x);
  int* nonzero_data = INTEGER(n_nonzero);
  int* max_feature_data = INTEGER(max_feature);
  double* umi_data = REAL(n_umi);
  double* mito_data = REAL(p_mito);
  double* loaded = REAL(y);
  double* fraction_data = REAL(frac_umis_max_feature);
  {
    std::vector<int> owned_mt(n_mt);
    for (R_xlen_t row = 0; row < n_mt; ++row) owned_mt[row] = mt_rows[row] - 1;
    for (int k = 0; k < n_cells; k ++) {
      int curr_max_feature = 0, n_nonzero_count = 0;
      double umi_count_curr_max_feature = 0, n_umi_count = 0;
      for (int row = 0; row < n_genes; ++row) loaded[row] = 0;
      for (int entry = pointers[k]; entry < pointers[k + 1]; ++entry)
        loaded[indices[entry]] = values[entry];
      for (int r = 0; r < n_genes; r ++) {
        n_umi_count += loaded[r];
        if (loaded[r] > 0.5) n_nonzero_count ++;
      }

      if (compute_p_mito) {
        double n_umi_count_mito = 0;
        for (int row : owned_mt) n_umi_count_mito += loaded[row];
        mito_data[k] = n_umi_count_mito/n_umi_count;
      }

      if (compute_max_feature) {
        for (int r = 0; r < n_genes; r ++) {
          if (loaded[r] > umi_count_curr_max_feature) {
            curr_max_feature = r;
            umi_count_curr_max_feature = loaded[r];
          }
        }
        max_feature_data[k] = curr_max_feature;
        fraction_data[k] = umi_count_curr_max_feature/n_umi_count;
      }

      nonzero_data[k] = n_nonzero_count;
      umi_data[k] = n_umi_count;
    }
  }
  return out;
}


// Count guide memberships per cell from a list of R one-based cell positions
// in the caller's n_cells view. Order/duplicates are not normalized: repeated
// memberships contribute repeatedly. Inputs have already been validated by the
// assignment path. Returns a new zero-initialized integer vector, no RNG/mutation.
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


// Threshold the STORED values of one CSR feature row (R one-based row_idx).
// Return R one-based source-cell positions in segment order for x>=threshold;
// explicit zeros can qualify at threshold<=0, but implicit zeros are not visited.
// Compare numeric fractions directly; no integer conversion or deduplication.
// Two passes allocate an exact result; source j/p/x and caller vectors unchanged.
// [[Rcpp::export]]
IntegerVector threshold_count_matrix(IntegerVector j, IntegerVector p, NumericVector x, int row_idx, double threshold) {
  IntegerVector out;
  int n_gte_threshold = 0, counter = 0, start = p[row_idx - 1], end = p[row_idx];
  for (int k = start; k < end; k ++) if (x[k] >= threshold) n_gte_threshold ++;
  out = IntegerVector(n_gte_threshold);
  for (int k = start; k < end; k ++) if (x[k] >= threshold) out[counter ++] = j[k] + 1;
  return out;
}
