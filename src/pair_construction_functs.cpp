#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
std::vector<int> compute_genes_within_distance(int midpoint, IntegerVector gene_tss_posits, int distance_threshold) {
  // the gene positions need not be sorted, so check every position
  std::vector<int> out;
  for (int i = 0; i < gene_tss_posits.size(); i ++) {
    if (abs(gene_tss_posits[i] - midpoint) <= distance_threshold) out.push_back(i + 1);
  }
  return(out);
}
