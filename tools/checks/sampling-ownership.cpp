// Developer-only instrumentation. The runner includes the actual candidate
// sampler source by absolute path and records both source and probe digests.
// Interpose only this translation unit's R allocation/preservation/finalizer
// calls; the installed package and frozen native oracle are untouched.
#define R_MakeExternalPtr refactor_controlled_make_external_ptr
#define R_GetCCallable refactor_controlled_get_callable
#define R_RegisterCFinalizerEx refactor_controlled_register_finalizer
#include <Rcpp.h>
#include <boost/random.hpp>
#include <cstring>
#undef R_MakeExternalPtr
#undef R_GetCCallable
#undef R_RegisterCFinalizerEx
#include SCEPTRE_SAMPLING_HEADER

// [[Rcpp::depends(BH)]]

extern "C" {
SEXP R_MakeExternalPtr(void*, SEXP, SEXP);
DL_FUNC R_GetCCallable(const char*, const char*);
void R_RegisterCFinalizerEx(SEXP, R_CFinalizer_t, Rboolean);
}

namespace {
int requested_phase = 0;
int observed_phase = 0;
bool address_was_null = false;

void fail_at_pointer_boundary(SEXP pointer, int phase) {
  if (requested_phase != phase || TYPEOF(pointer) != EXTPTRSXP) return;
  observed_phase = phase;
  address_was_null = R_ExternalPtrAddr(pointer) == nullptr;
  requested_phase = 0;
  // On an old bank-first implementation, use a C++ exception to permit RAII
  // cleanup instead of deliberately leaking its live native allocation.
  if (!address_was_null) Rcpp::stop("Bank was allocated before the R boundary.");
  Rf_error("Controlled R pointer-boundary failure.");
}

SEXP controlled_precious_preserve(SEXP object) {
  fail_at_pointer_boundary(object, 2);
  typedef SEXP (*Preserve)(SEXP);
  static Preserve preserve = reinterpret_cast<Preserve>(
    R_GetCCallable("Rcpp", "Rcpp_precious_preserve"));
  return preserve(object);
}
}

extern "C" SEXP refactor_controlled_make_external_ptr(void* address, SEXP tag, SEXP protection) {
  if (requested_phase == 1) {
    observed_phase = 1;
    address_was_null = address == nullptr;
    requested_phase = 0;
    if (!address_was_null) Rcpp::stop("Bank was allocated before the R boundary.");
    Rf_error("Controlled R pointer-boundary failure.");
  }
  return R_MakeExternalPtr(address, tag, protection);
}

extern "C" DL_FUNC refactor_controlled_get_callable(const char* package, const char* name) {
  if (std::strcmp(package, "Rcpp") == 0 && std::strcmp(name, "Rcpp_precious_preserve") == 0)
    return reinterpret_cast<DL_FUNC>(&controlled_precious_preserve);
  return R_GetCCallable(package, name);
}

extern "C" void refactor_controlled_register_finalizer(SEXP pointer, R_CFinalizer_t finalizer, Rboolean on_exit) {
  fail_at_pointer_boundary(pointer, 3);
  R_RegisterCFinalizerEx(pointer, finalizer, on_exit);
}

namespace actual_candidate_samplers {
#include SCEPTRE_SAMPLING_SOURCE
}

// [[Rcpp::export]]
void refactor_sampling_pointer_failure(int sampler, int phase) {
  requested_phase = phase;
  observed_phase = 0;
  address_was_null = false;
  if (sampler == 0) actual_candidate_samplers::fisher_yates_sampler(3, 1, 2);
  else if (sampler == 1) actual_candidate_samplers::hybrid_fisher_iwor_sampler(2, 1, 2, 2);
  else if (sampler == 2) actual_candidate_samplers::sample_crt_indices(Rcpp::NumericVector::create(0.5), 2);
  else Rcpp::stop("Unknown sampler probe.");
  requested_phase = 0;
}

// [[Rcpp::export]]
Rcpp::List refactor_sampling_failure_state() {
  return Rcpp::List::create(Rcpp::Named("phase") = observed_phase,
                           Rcpp::Named("address_was_null") = address_was_null);
}
