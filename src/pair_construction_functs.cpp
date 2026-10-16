#include <Rcpp.h>
#include <cstdint>
#include <limits>
using namespace Rcpp;


// Supplied chromosome positions need not be sorted. Return their one-based
// input positions satisfying the inclusive distance predicate, in input order.
// Widen before subtraction so valid signed integer coordinates cannot overflow.
// Negative thresholds retain an empty result; missing/unrepresentable native
// coordinates fail clearly. Inputs are read only and no RNG is used.
// Count before allocating the exact R result: no native-owned result buffer is
// alive across R allocation/return wrapping, and no full-length scratch is needed.
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
