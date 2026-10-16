#ifndef SCEPTRE_SAMPLING_SAFETY_H
#define SCEPTRE_SAMPLING_SAFETY_H

#include <Rcpp.h>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

// Preconditions for existing int-indexed samplers. Check before arithmetic or
// allocation; these helpers do not choose smaller budgets or change draw order.
namespace sceptre {
namespace sampling {
inline void nonnegative(int value, const char* name) {
  if (value < 0) Rcpp::stop("%s must be nonnegative and not NA.", name);
}

inline int sum(int left, int right, const char* name) {
  nonnegative(left, name);
  nonnegative(right, name);
  const int64_t total = static_cast<int64_t>(left) + right;
  if (total > std::numeric_limits<int>::max())
    Rcpp::stop("%s exceeds native integer capacity.", name);
  return static_cast<int>(total);
}

inline void matrix_extent(int rows, int columns) {
  nonnegative(rows, "Matrix rows");
  nonnegative(columns, "Matrix columns");
  const uint64_t size = static_cast<uint64_t>(rows) * columns;
  if (size > static_cast<uint64_t>(R_XLEN_T_MAX) ||
      size > std::vector<int>().max_size())
    Rcpp::stop("Requested matrix exceeds allocation capacity.");
}

inline void bank_extent(int rows, int width) {
  nonnegative(rows, "B");
  nonnegative(width, "Draw width");
  if (static_cast<size_t>(rows) > std::vector<std::vector<int>>().max_size() ||
      static_cast<size_t>(width) > std::vector<int>().max_size() ||
      static_cast<uint64_t>(rows) * width >
        std::numeric_limits<size_t>::max() / sizeof(int))
    Rcpp::stop("Requested draw bank exceeds allocation capacity.");
}

// The finalizing NULL-address pointer must already be protected before the
// native bank is allocated. Creating/preserving/registering an XPtr may raise
// an R error that bypasses C++ destructors; this final transfer cannot allocate.
inline SEXP register_bank(Rcpp::XPtr<std::vector<std::vector<int>>>& pointer,
                          std::unique_ptr<std::vector<std::vector<int>>>& bank) {
  R_SetExternalPtrAddr(pointer, bank.get());
  bank.release();
  return pointer;
}
} // namespace sampling
} // namespace sceptre
#endif
