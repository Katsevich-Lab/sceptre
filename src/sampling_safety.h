#ifndef SCEPTRE_SAMPLING_SAFETY_H
#define SCEPTRE_SAMPLING_SAFETY_H

#include <Rcpp.h>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

// Check sampler sizes before arithmetic or allocation. Invalid requests error;
// these helpers never reduce draw counts or alter the random sampling order.
namespace sceptre {
namespace sampling {
// Reject negative counts, including R's negative NA_INTEGER marker, before using
// them in sizes or arithmetic. Zero is allowed where the sampler supports it.
inline void nonnegative(int value, const char* name) {
  if (value < 0) Rcpp::stop("%s must be nonnegative and not NA.", name);
}

// Add two nonnegative counts using a wider integer, then check that the sum fits
// an int. Use for stage offsets and total cell counts; error rather than cap it.
inline int sum(int left, int right, const char* name) {
  nonnegative(left, name);
  nonnegative(right, name);
  const int64_t total = static_cast<int64_t>(left) + right;
  if (total > std::numeric_limits<int>::max())
    Rcpp::stop("%s exceeds native integer capacity.", name);
  return static_cast<int>(total);
}

// Check that rows * columns fits R's length limit and C++ integer-vector size
// limits before allocating the R matrix. This checks representable size,
// not whether the machine has enough free memory.
inline void matrix_extent(int rows, int columns) {
  nonnegative(rows, "Matrix rows");
  nonnegative(columns, "Matrix columns");
  const uint64_t size = static_cast<uint64_t>(rows) * columns;
  if (size > static_cast<uint64_t>(R_XLEN_T_MAX) ||
      size > std::vector<int>().max_size())
    Rcpp::stop("Requested matrix exceeds allocation capacity.");
}

// Check that the number of draws, entries per draw and total integer bytes fit
// the C++ container and size limits. Empty draws are allowed. No memory for
// draws is allocated here; passing this check does not guarantee enough RAM.
inline void bank_extent(int rows, int width) {
  nonnegative(rows, "B");
  nonnegative(width, "Draw width");
  if (static_cast<size_t>(rows) > std::vector<std::vector<int>>().max_size() ||
      static_cast<size_t>(width) > std::vector<int>().max_size() ||
      static_cast<uint64_t>(rows) * width >
        std::numeric_limits<size_t>::max() / sizeof(int))
    Rcpp::stop("Requested draw bank exceeds allocation capacity.");
}

// Attach completed draws to an R external pointer that is already protected and
// has its cleanup handler registered. Pointer creation can raise an R error,
// which skips C++ destructors, so it must happen before allocating the draws.
// This final ownership transfer must not allocate or call an error-producing R
// operation; otherwise the draws could be leaked between their two owners.
inline SEXP register_bank(Rcpp::XPtr<std::vector<std::vector<int>>>& pointer,
                          std::unique_ptr<std::vector<std::vector<int>>>& bank) {
  R_SetExternalPtrAddr(pointer, bank.get());
  bank.release();
  return pointer;
}
} // namespace sampling
} // namespace sceptre
#endif
