#ifndef SCEPTRE_RESULT_LIFETIME_H
#define SCEPTRE_RESULT_LIFETIME_H

#include <Rcpp.h>

namespace sceptre {
// Create and protect R's error-continuation token before allocating C++ results.
// R errors can skip C++ destructors; the token lets a later R_UnwindProtect call
// run cleanup even on that path. If ordinary C++ computation throws, this guard
// removes the protection normally. release() transfers that responsibility to
// the caller without calling R. End the guard's scope before R_UnwindProtect,
// and UNPROTECT once after successful conversion. On an R error, R restores
// its protection stack itself.
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
