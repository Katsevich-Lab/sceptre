#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
std::vector<int> compute_genes_within_distance(int midpoint, IntegerVector gene_tss_posits, int distance_threshold) {
  std::vector<int> out;
  bool within_distance_threshold = false, hit_distance_threshold = false;
  for (int i = 0; i < gene_tss_posits.size(); i ++) {
    within_distance_threshold = abs(gene_tss_posits[i] - midpoint) <= distance_threshold;
    if (within_distance_threshold) {
      hit_distance_threshold = true;
      out.push_back(i + 1);
    }
    if (hit_distance_threshold && !within_distance_threshold) break;
  }
  return(out);
}
