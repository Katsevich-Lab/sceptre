#include <Rcpp.h>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace Rcpp;

namespace {
// Trivial views use R's transient arena, so an R error during preparation never
// strands a native owner. Protected normalized R inputs retain the viewed data.
struct PositionView {
  const int* data;
  int size;
};

PositionView prepare_positions(SEXP positions, int limit, const char* name) {
  const R_xlen_t size = Rf_xlength(positions);
  if (size > std::numeric_limits<int>::max())
    stop("%s exceeds native coordinate capacity.", name);
  const int* data = INTEGER(positions); // Materialize ALTREP before native ownership.
  for (R_xlen_t index = 0; index < size; ++index)
    if (data[index] < 1 || data[index] > limit)
      stop("%s contains a position outside its one-based coordinate view.", name);
  return {data, static_cast<int>(size)};
}

PositionView* prepare_groups(const List& groups, int n_cells, SEXP normalized) {
  PositionView* views = reinterpret_cast<PositionView*>(R_alloc(groups.size(), sizeof(PositionView)));
  for (R_xlen_t group = 0; group < groups.size(); ++group) {
    Shield<SEXP> positions(Rf_coerceVector(VECTOR_ELT(groups, group), INTSXP));
    views[group] = prepare_positions(positions, n_cells, "Analysis group");
    SET_VECTOR_ELT(normalized, group, positions);
  }
  return views;
}

int checked_group_count(const List& groups) {
  if (groups.size() > std::numeric_limits<int>::max())
    stop("Group count exceeds native coordinate capacity.");
  return static_cast<int>(groups.size());
}

int response_count(const IntegerVector& p, int n_cells_orig,
                   const IntegerVector& cells) {
  if (n_cells_orig < 0 || p.size() < 1 ||
      p.size() - 1 > std::numeric_limits<int>::max() ||
      cells.size() > std::numeric_limits<int>::max())
    stop("Sparse count dimensions exceed their native coordinate capacity.");
  return static_cast<int>(p.size() - 1);
}

// Pure C++ helpers make no R calls; native allocation exceptions unwind their
// owners. Repetitions/order survive the sole one-based-to-zero-based copy.
std::vector<int> copy_zero_based(PositionView positions) {
  std::vector<int> owned(positions.size);
  for (int index = 0; index < positions.size; ++index)
    owned[index] = positions.data[index] - 1;
  return owned;
}

std::vector<std::vector<int>> copy_zero_based_groups(const PositionView* groups,
                                                    int n_groups) {
  std::vector<std::vector<int>> owned(n_groups);
  for (int group = 0; group < n_groups; ++group)
    owned[group] = copy_zero_based(groups[group]);
  return owned;
}

// CSR slots are valid zero-based Matrix coordinates. Presence includes stored
// zeros/fractions. Source and analysis scratch belong only to this native scope.
void load_stored_positions(const int* j, const int* p, int response,
                           std::vector<bool>& source,
                           std::vector<bool>& analysis,
                           const std::vector<int>& cells) {
  std::fill(source.begin(), source.end(), false);
  for (int entry = p[response]; entry < p[response + 1]; ++entry)
    source[j[entry]] = true;
  for (size_t cell = 0; cell < cells.size(); ++cell)
    analysis[cell] = source[cells[cell]];
}
}

// Fused all-response summaries and requested pair counts from valid CSR slots.
// R one-based groups index the ordered cells_in_use analysis view directly;
// requests are response-sorted and storage restores their original row order.
// NT totals SUM groups, including overlap; complement totals count analysis
// cells. Coercion/materialization and named R results precede owned native
// copies. Their compute-only scope ends before return, including on exceptions.
// No input mutation, count-value loading, eligibility, RNG, or file effects.
// [[Rcpp::export]]
List sparse_pair_counts(IntegerVector j, IntegerVector p, int n_cells_orig,
                        List grna_group_idxs, List nt_groups_analysis,
                        IntegerVector to_analyze_response_idxs,
                        IntegerVector to_analyze_grna_idxs,
                        bool control_group_complement,
                        IntegerVector cells_in_use) {
  const int n_responses = response_count(p, n_cells_orig, cells_in_use);
  const int n_cells = static_cast<int>(cells_in_use.size());
  if (to_analyze_response_idxs.size() != to_analyze_grna_idxs.size() ||
      to_analyze_response_idxs.size() > std::numeric_limits<int>::max())
    stop("Requested response/group positions must have matching native lengths.");
  PositionView cells = prepare_positions(cells_in_use, n_cells_orig, "Source cells");
  const int n_groups = checked_group_count(grna_group_idxs);
  const int n_nt_groups = checked_group_count(nt_groups_analysis);
  Shield<SEXP> normalized_groups(Rf_allocVector(VECSXP, n_groups));
  Shield<SEXP> normalized_nt(Rf_allocVector(VECSXP, n_nt_groups));
  PositionView* group_views = prepare_groups(grna_group_idxs, n_cells, normalized_groups);
  PositionView* nt_views = prepare_groups(nt_groups_analysis, n_cells, normalized_nt);
  PositionView responses = prepare_positions(to_analyze_response_idxs, n_responses,
                                             "Requested response");
  PositionView pairs = prepare_positions(to_analyze_grna_idxs, n_groups, "Requested group");
  if (responses.size > 1 && !std::is_sorted(responses.data, responses.data + responses.size))
    stop("Sparse pair requests must be in response scan order.");
  const int* columns = INTEGER(j);
  const int* pointers = INTEGER(p);

  IntegerVector totals(n_responses), treatment_counts(responses.size), control_counts(responses.size);
  IntegerMatrix nt_counts(n_nt_groups, n_responses);
  List result = List::create(Named("n_nonzero_mat") = nt_counts, Named("n_nonzero_tot") = totals,
                            Named("n_nonzero_trt") = treatment_counts, Named("n_nonzero_cntrl") = control_counts);
  int* total_data = INTEGER(totals);
  int* treatment_data = INTEGER(treatment_counts);
  int* control_data = INTEGER(control_counts);
  int* nt_data = INTEGER(nt_counts);
  {
    const auto owned_cells = copy_zero_based(cells);
    const auto groups = copy_zero_based_groups(group_views, n_groups);
    const auto nt_groups = copy_zero_based_groups(nt_views, n_nt_groups);
    const auto requested_responses = copy_zero_based(responses);
    const auto requested_groups = copy_zero_based(pairs);
    std::vector<bool> stored_in_source(n_cells_orig), stored_in_analysis(n_cells);
    int pair_row = 0;
    for (int response = 0; response < n_responses; ++response) {
      load_stored_positions(columns, pointers, response, stored_in_source, stored_in_analysis, owned_cells);
      int64_t nt_total = 0;
      for (int group = 0; group < n_nt_groups; ++group) {
        int count = 0;
        for (int cell : nt_groups[group]) if (stored_in_analysis[cell]) ++count;
        nt_data[group + static_cast<R_xlen_t>(response) * n_nt_groups] = count;
        nt_total += count;
      }
      if (control_group_complement) {
        for (bool stored : stored_in_analysis) if (stored) ++total_data[response];
      } else {
        if (nt_total > std::numeric_limits<int>::max())
          throw std::overflow_error("NT group sum exceeds integer count capacity.");
        total_data[response] = static_cast<int>(nt_total);
      }
      while (pair_row < responses.size && requested_responses[pair_row] == response) {
        int treated = 0;
        for (int cell : groups[requested_groups[pair_row]])
          if (stored_in_analysis[cell]) ++treated;
        treatment_data[pair_row] = treated;
        control_data[pair_row] = control_group_complement ? total_data[response] - treated : total_data[response];
        ++pair_row;
      }
    }
  }
  return result;
}

// Group x all-response stored-position counts, using the same private loader.
// Normalized inputs and the complete matrix exist before owned native copies;
// the compute-only scope ends before return. Empty views retain zero dimensions
// and counts; storage selects/repeats columns later. No input mutation,
// count-value loading, RNG, or file effects.
// [[Rcpp::export]]
IntegerMatrix sparse_group_counts(IntegerVector j, IntegerVector p, int n_cells_orig,
                                  List groups, IntegerVector cells_in_use) {
  const int n_responses = response_count(p, n_cells_orig, cells_in_use);
  const int n_cells = static_cast<int>(cells_in_use.size());
  PositionView cells = prepare_positions(cells_in_use, n_cells_orig, "Source cells");
  const int n_groups = checked_group_count(groups);
  Shield<SEXP> normalized_groups(Rf_allocVector(VECSXP, n_groups));
  PositionView* group_views = prepare_groups(groups, n_cells, normalized_groups);
  const int* columns = INTEGER(j);
  const int* pointers = INTEGER(p);
  IntegerMatrix counts(n_groups, n_responses);
  int* count_data = INTEGER(counts);
  {
    const auto owned_cells = copy_zero_based(cells);
    const auto owned_groups = copy_zero_based_groups(group_views, n_groups);
    std::vector<bool> stored_in_source(n_cells_orig), stored_in_analysis(n_cells);
    for (int response = 0; response < n_responses; ++response) {
      load_stored_positions(columns, pointers, response, stored_in_source, stored_in_analysis, owned_cells);
      for (int group = 0; group < n_groups; ++group) {
        int count = 0;
        for (int cell : owned_groups[group]) if (stored_in_analysis[cell]) ++count;
        count_data[group + static_cast<R_xlen_t>(response) * n_groups] = count;
      }
    }
  }
  return counts;
}
