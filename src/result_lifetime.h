#ifndef SCEPTRE_RESULT_LIFETIME_H
#define SCEPTRE_RESULT_LIFETIME_H

#include <Rcpp.h>

namespace sceptre {
// Protect the R continuation before native result ownership begins. A C++
// compute failure unwinds this guard normally; release transfers only its
// protection obligation, without an R call. End this guard's scope before
// entering R_UnwindProtect, then UNPROTECT once after the protected leaf returns.
// The leaf's raw-R-error path restores the R protection stack itself.
class ResultContinuation {
 public:
  ResultContinuation() : token_(PROTECT(R_MakeUnwindCont())), active_(true) {}
  ~ResultContinuation() { if (active_) UNPROTECT(1); }
  SEXP release() noexcept { active_ = false; return token_; }
 private:
  SEXP token_;
  bool active_;
  ResultContinuation(const ResultContinuation&) = delete;
  ResultContinuation& operator=(const ResultContinuation&) = delete;
};
} // namespace sceptre

#endif
