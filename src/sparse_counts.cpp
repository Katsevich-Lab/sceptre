#include <Rcpp.h>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace Rcpp;

namespace {
// Keep these small position records in R-managed temporary memory. If R errors
// while preparing inputs, R can reclaim them; there is no C++ buffer left to
// leak. Keep the R vectors protected while these records point into them.
struct PositionView {
  const int* data;
  int size;
};

// Read an integer vector of positions starting at 1. Check that positions are
// between 1 and limit and that the length fits an int. Keep order and repeats;
// do not subtract from the input. The caller keeps this R vector protected
// until C++ has finished using its data pointer.
PositionView prepare_positions(SEXP positions, int limit, const char* name) {
  const R_xlen_t size = Rf_xlength(positions);
  if (size > std::numeric_limits<int>::max())
    stop("%s exceeds native coordinate capacity.", name);
  const int* data = INTEGER(positions); // Force any lazy R data to load before creating C++ buffers.
  for (R_xlen_t index = 0; index < size; ++index)
    if (data[index] < 1 || data[index] > limit)
      stop("%s contains a position outside its one-based coordinate view.", name);
  return {data, static_cast<int>(size)};
}

// Convert each group to integers and keep it in the protected normalized R list
// before allocating C++ buffers. An R error during conversion can skip C++
// destructors; use R-managed temporary memory for these position records.
// Group members are positions among the selected analysis cells.
PositionView* prepare_groups(const List& groups, int n_cells, SEXP normalized) {
  PositionView* views = reinterpret_cast<PositionView*>(R_alloc(groups.size(), sizeof(PositionView)));
  for (R_xlen_t group = 0; group < groups.size(); ++group) {
    Shield<SEXP> positions(Rf_coerceVector(VECTOR_ELT(groups, group), INTSXP));
    views[group] = prepare_positions(positions, n_cells, "Analysis group");
    SET_VECTOR_ELT(normalized, group, positions);
  }
  return views;
}

// Check that the group count fits an int before using it as a matrix dimension.
// No groups is valid and gives an output with zero rows.
int checked_group_count(const List& groups) {
  if (groups.size() > std::numeric_limits<int>::max())
    stop("Group count exceeds native coordinate capacity.");
  return static_cast<int>(groups.size());
}

// Get the number of gene rows from the row-pointer vector and check that gene
// and cell counts fit integer positions. The storage caller must already supply
// valid compressed sparse rows (CSR); this is not a full sparse-format check.
int response_count(const IntegerVector& p, int n_cells_orig,
                   const IntegerVector& cells) {
  if (n_cells_orig < 0 || p.size() < 1 ||
      p.size() - 1 > std::numeric_limits<int>::max() ||
      cells.size() > std::numeric_limits<int>::max())
    stop("Sparse count dimensions exceed their native coordinate capacity.");
  return static_cast<int>(p.size() - 1);
}

// The helpers below copy and convert positions using C++ only, without R calls.
// C++ exceptions free their buffers normally. Subtract one during the copy,
// keeping membership order and repeated cells.
std::vector<int> copy_zero_based(PositionView positions) {
  std::vector<int> owned(positions.size);
  for (int index = 0; index < positions.size; ++index)
    owned[index] = positions.data[index] - 1;
  return owned;
}

// After R input preparation, make private C++ copies of all groups and convert
// positions from starting at 1 to starting at 0. Keep empty groups and repeated
// members, including repeats contributing to NT totals. Only C++ operations may
// run while these buffers need C++ destructors for cleanup.
std::vector<std::vector<int>> copy_zero_based_groups(const PositionView* groups,
                                                    int n_groups) {
  std::vector<std::vector<int>> owned(n_groups);
  for (int group = 0; group < n_groups; ++group)
    owned[group] = copy_zero_based(groups[group]);
  return owned;
}

// Read valid compressed sparse rows using Matrix's zero-based row/column data.
// Count stored positions, including stored zeros and fractions, without loading
// their values. Temporary cell mappings belong only to this C++ calculation.
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

// Count requested gene-guide pairs and collect NT summaries in the same scan.
// Groups use positions starting at 1 within the ordered cells_in_use list.
// Requests arrive sorted by gene; R storage code restores the requested order.
// NT totals add group counts, counting overlaps again; complement totals count
// selected cells. These are counts of stored entries, not positive values.
// Prepare R inputs and result objects before creating private C++ position copies.
// End the C++ calculation scope before returning to R, freeing copies on success
// or C++ error. Leave inputs unchanged and apply QC thresholds separately.
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

// Return stored-entry counts with guide groups in rows and all genes in columns.
// Use the same cell-mapping helper as pair counts. Empty selections retain the
// appropriate zero dimensions/counts; R selects and repeats gene columns later.
// Prepare R inputs and the result matrix first, then create C++ copies and finish
// using them before returning. Caller inputs are unchanged; count values need
// not be loaded because this operation counts stored positions.
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
