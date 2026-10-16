#include <Rcpp.h>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

using namespace Rcpp;

namespace {
// CSR slots are valid zero-based Matrix coordinates. Count presence of stored
// positions, including stored zeros/fractions, rather than loaded-value > 0.
// The source/analysis scratch vectors belong to this call; cells is the owned
// zero-based analysis-to-source map and may retain requested repetitions.
void load_stored_positions(const IntegerVector& j, const IntegerVector& p,
                           int response, std::vector<bool>& source,
                           std::vector<bool>& analysis,
                           const std::vector<int>& cells) {
  std::fill(source.begin(), source.end(), false);
  for (int entry = p[response]; entry < p[response + 1]; ++entry)
    source[j[entry]] = true;
  for (size_t cell = 0; cell < cells.size(); ++cell)
    analysis[cell] = source[cells[cell]];
}

// Rcpp copied these position vectors on entry. Convert exactly once without
// sorting/deduplication: repetitions and each group's order remain meaningful.
void zero_base(std::vector<int>& positions, int limit, const char* name) {
  for (int& position : positions) {
    if (position < 1 || position > limit)
      stop("%s contains a position outside its one-based coordinate view.", name);
    --position;
  }
}

int response_count(const IntegerVector& p, int n_cells_orig,
                   const std::vector<int>& cells) {
  if (n_cells_orig < 0 || p.size() < 1 ||
      p.size() - 1 > std::numeric_limits<int>::max() ||
      cells.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
    stop("Sparse count dimensions exceed their native coordinate capacity.");
  return static_cast<int>(p.size() - 1);
}

int zero_base_groups(std::vector<std::vector<int>>& groups, int n_cells) {
  if (groups.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
    stop("Group count exceeds native coordinate capacity.");
  for (auto& group : groups) {
    if (group.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
      stop("Group membership exceeds integer count capacity.");
    zero_base(group, n_cells, "Analysis group");
  }
  return static_cast<int>(groups.size());
}
}

// Fused all-response summaries and requested pair counts from valid CSR slots.
// All positions except j/p arrive R one-based. Both treatment and NT groups
// index the ordered cells_in_use analysis view directly; no NT-local round trip
// occurs. Requests must be response-sorted; storage restores the original rows.
// Owned coordinate conversion leaves every R input unchanged. Returns integer
// NT-group x all-response counts, all-response totals, and requested pair counts.
// Complement totals count analysis cells; NT totals deliberately SUM NT groups,
// including overlapping memberships. No eligibility rule, RNG or file effects.
// [[Rcpp::export]]
List sparse_pair_counts(IntegerVector j, IntegerVector p, int n_cells_orig,
                        std::vector<std::vector<int>> grna_group_idxs,
                        std::vector<std::vector<int>> nt_groups_analysis,
                        std::vector<int> to_analyze_response_idxs,
                        std::vector<int> to_analyze_grna_idxs,
                        bool control_group_complement,
                        std::vector<int> cells_in_use) {
  const int n_responses = response_count(p, n_cells_orig, cells_in_use);
  const int n_cells = cells_in_use.size();
  if (to_analyze_response_idxs.size() != to_analyze_grna_idxs.size() ||
      to_analyze_response_idxs.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
    stop("Requested response/group positions must have matching native lengths.");
  zero_base(cells_in_use, n_cells_orig, "Source cells");
  const int n_groups = zero_base_groups(grna_group_idxs, n_cells);
  const int n_nt_groups = zero_base_groups(nt_groups_analysis, n_cells);
  zero_base(to_analyze_response_idxs, n_responses, "Requested response");
  zero_base(to_analyze_grna_idxs, n_groups, "Requested group");
  if (!std::is_sorted(to_analyze_response_idxs.begin(), to_analyze_response_idxs.end()))
    stop("Sparse pair requests must be in response scan order.");

  const int n_pairs = to_analyze_response_idxs.size();
  IntegerVector totals(n_responses), treatment_counts(n_pairs), control_counts(n_pairs);
  IntegerMatrix nt_counts(n_nt_groups, n_responses);
  std::vector<bool> stored_in_source(n_cells_orig), stored_in_analysis(n_cells);
  int pair_row = 0;

  for (int response = 0; response < n_responses; ++response) {
    load_stored_positions(j, p, response, stored_in_source, stored_in_analysis, cells_in_use);
    int64_t nt_total = 0;
    for (int group = 0; group < n_nt_groups; ++group) {
      int count = 0;
      for (int cell : nt_groups_analysis[group]) if (stored_in_analysis[cell]) ++count;
      nt_counts(group, response) = count;
      nt_total += count;
    }
    if (control_group_complement) {
      for (bool stored : stored_in_analysis) if (stored) ++totals[response];
    } else {
      if (nt_total > std::numeric_limits<int>::max()) stop("NT group sum exceeds integer count capacity.");
      totals[response] = static_cast<int>(nt_total);
    }
    while (pair_row < n_pairs && to_analyze_response_idxs[pair_row] == response) {
      int treated = 0;
      for (int cell : grna_group_idxs[to_analyze_grna_idxs[pair_row]])
        if (stored_in_analysis[cell]) ++treated;
      treatment_counts[pair_row] = treated;
      control_counts[pair_row] = control_group_complement ? totals[response] - treated : totals[response];
      ++pair_row;
    }
  }
  return List::create(Named("n_nonzero_mat") = nt_counts, Named("n_nonzero_tot") = totals,
                      Named("n_nonzero_trt") = treatment_counts, Named("n_nonzero_cntrl") = control_counts);
}

// Group x all-response stored-position counts, using the same private loader as
// pair summaries. Rcpp owns group/cell vectors; R one-based analysis positions
// become zero-based once. Empty groups/views retain zero counts and dimensions;
// the R storage adapter selects/repeats requested response columns afterwards.
// No caller mutation, count-value loading, RNG, or file effects.
// [[Rcpp::export]]
IntegerMatrix sparse_group_counts(IntegerVector j, IntegerVector p, int n_cells_orig,
                                   std::vector<std::vector<int>> groups,
                                   std::vector<int> cells_in_use) {
  const int n_responses = response_count(p, n_cells_orig, cells_in_use);
  const int n_cells = cells_in_use.size();
  zero_base(cells_in_use, n_cells_orig, "Source cells");
  const int n_groups = zero_base_groups(groups, n_cells);
  IntegerMatrix counts(n_groups, n_responses);
  std::vector<bool> stored_in_source(n_cells_orig), stored_in_analysis(n_cells);
  for (int response = 0; response < n_responses; ++response) {
    load_stored_positions(j, p, response, stored_in_source, stored_in_analysis, cells_in_use);
    for (size_t group = 0; group < groups.size(); ++group) {
      int count = 0;
      for (int cell : groups[group]) if (stored_in_analysis[cell]) ++count;
      counts(group, response) = count;
    }
  }
  return counts;
}
