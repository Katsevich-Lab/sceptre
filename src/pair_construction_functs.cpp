#include <Rcpp.h>
#include <cstdint>
#include <limits>
using namespace Rcpp;


// Find supplied genomic positions within the inclusive distance threshold,
// even if they are unsorted. Return their positions in the input, starting at 1,
// in input order. Use wider arithmetic before subtracting genomic coordinates
// to avoid integer overflow. Negative thresholds give no matches; missing or
// out-of-range integer coordinates error. Count matches before allocating the R
// output, so no C++ result buffer can leak if R allocation fails.
// [[Rcpp::export]]
IntegerVector compute_genes_within_distance(int midpoint, IntegerVector gene_tss_posits, int distance_threshold) {
  if (midpoint == NA_INTEGER || distance_threshold == NA_INTEGER ||
      gene_tss_posits.size() > std::numeric_limits<int>::max())
    stop("Cis midpoint, threshold and response positions must fit native integer coordinates.");
  const int n = static_cast<int>(gene_tss_posits.size());
  const int64_t distance = distance_threshold;
  const auto within_distance = [midpoint, distance](int position) {
    if (position == NA_INTEGER) stop("Cis response positions cannot be missing.");
    const int64_t delta = static_cast<int64_t>(position) - midpoint;
    return delta >= -distance && delta <= distance;
  };
  int n_matches = 0;
  for (int row = 0; row < n; ++row) if (within_distance(gene_tss_posits[row])) ++n_matches;
  IntegerVector result(n_matches);
  int next = 0;
  for (int row = 0; row < n; ++row)
    if (within_distance(gene_tss_posits[row])) result[next++] = row + 1;
  return result;
}
