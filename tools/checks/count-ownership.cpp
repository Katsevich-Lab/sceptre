// Test memory cleanup in the package's actual count functions. Wrap R allocation
// calls before including Rcpp templates, so the test can trigger errors during
// vector/matrix creation, named-list construction and Rcpp object protection.
#define Rf_allocVector refactor_count_alloc_vector
#define Rf_allocMatrix refactor_count_alloc_matrix
#define R_alloc refactor_count_r_alloc
#define Rf_coerceVector refactor_count_coerce
#define INTEGER refactor_count_integer
#define REAL refactor_count_real
#define Rf_mkChar refactor_count_mkchar
#define Rf_setAttrib refactor_count_setattr
#define R_GetCCallable refactor_count_get_callable
#include <Rcpp.h>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>
#undef Rf_allocVector
#undef Rf_allocMatrix
#undef R_alloc
#undef Rf_coerceVector
#undef INTEGER
#undef REAL
#undef Rf_mkChar
#undef Rf_setAttrib
#undef R_GetCCallable

extern "C" {
SEXP Rf_allocVector(SEXPTYPE, R_xlen_t);
SEXP Rf_allocMatrix(SEXPTYPE, int, int);
char* R_alloc(R_SIZE_T, int);
SEXP Rf_coerceVector(SEXP, SEXPTYPE);
int* INTEGER(SEXP);
double* REAL(SEXP);
SEXP Rf_mkChar(const char*);
SEXP Rf_setAttrib(SEXP, SEXP, SEXP);
DL_FUNC R_GetCCallable(const char*, const char*);
}

namespace {
bool observing = false;
int requested_boundary = 0, requested_native = 0, native_error_kind = 0;
int boundary_count = 0, native_count = 0, live_blocks = 0, max_live_blocks = 0;
int failed_boundary = 0, failed_native = 0, unsafe_boundaries = 0;
int observed_kinds[4096], observed_live[4096];

void r_boundary(int kind) {
  if (!observing) return;
  if (boundary_count == 4096) throw std::runtime_error("Probe boundary capacity exceeded.");
  observed_kinds[boundary_count] = kind;
  observed_live[boundary_count++] = live_blocks;
  if (live_blocks != 0) {
    ++unsafe_boundaries;
    // Avoid deliberately leaking the defective path while making it fail.
    throw std::runtime_error("R boundary observed with live native count ownership.");
  }
  if (boundary_count == requested_boundary) {
    failed_boundary = boundary_count;
    requested_boundary = 0;
    Rf_error("Controlled count R-boundary failure.");
  }
}

SEXP count_precious_preserve(SEXP object) {
  r_boundary(9);
  typedef SEXP (*Preserve)(SEXP);
  static Preserve preserve = reinterpret_cast<Preserve>(R_GetCCallable("Rcpp", "Rcpp_precious_preserve"));
  return preserve(object);
}
}

extern "C" SEXP refactor_count_alloc_vector(SEXPTYPE type, R_xlen_t length) {
  r_boundary(1); return Rf_allocVector(type, length);
}
extern "C" SEXP refactor_count_alloc_matrix(SEXPTYPE type, int rows, int columns) {
  r_boundary(2); return Rf_allocMatrix(type, rows, columns);
}
extern "C" char* refactor_count_r_alloc(R_SIZE_T count, int width) {
  r_boundary(3); return R_alloc(count, width);
}
extern "C" SEXP refactor_count_coerce(SEXP object, SEXPTYPE type) {
  r_boundary(4); return Rf_coerceVector(object, type);
}
extern "C" int* refactor_count_integer(SEXP object) {
  r_boundary(5); return INTEGER(object);
}
extern "C" double* refactor_count_real(SEXP object) {
  r_boundary(6); return REAL(object);
}
extern "C" SEXP refactor_count_mkchar(const char* string) {
  r_boundary(7); return Rf_mkChar(string);
}
extern "C" SEXP refactor_count_setattr(SEXP object, SEXP name, SEXP value) {
  r_boundary(8);
  if (name == R_DimSymbol) r_boundary(11);
  return Rf_setAttrib(object, name, value);
}
extern "C" DL_FUNC refactor_count_get_callable(const char* package, const char* name) {
  r_boundary(10);
  if (std::strcmp(package, "Rcpp") == 0 && std::strcmp(name, "Rcpp_precious_preserve") == 0)
    return reinterpret_cast<DL_FUNC>(&count_precious_preserve);
  return R_GetCCallable(package, name);
}

template<class T> struct RefactorCountAllocator {
  using value_type = T;
  RefactorCountAllocator() noexcept = default;
  template<class U> RefactorCountAllocator(const RefactorCountAllocator<U>&) noexcept {}
  T* allocate(std::size_t count) {
    ++native_count;
    if (observing && native_count == requested_native) {
      failed_native = native_count;
      if (native_error_kind == 1) throw std::overflow_error("Controlled native count overflow.");
      throw std::bad_alloc();
    }
    T* pointer = std::allocator<T>{}.allocate(count);
    ++live_blocks;
    max_live_blocks = std::max(max_live_blocks, live_blocks);
    return pointer;
  }
  void deallocate(T* pointer, std::size_t count) noexcept {
    --live_blocks;
    std::allocator<T>{}.deallocate(pointer, count);
  }
  template<class U> bool operator==(const RefactorCountAllocator<U>&) const noexcept { return true; }
  template<class U> bool operator!=(const RefactorCountAllocator<U>&) const noexcept { return false; }
};
template<class T> using refactor_count_vector = std::vector<T, RefactorCountAllocator<T>>;

int refactor_count_overflow_limit = std::numeric_limits<int>::max();
#define Rf_allocVector refactor_count_alloc_vector
#define Rf_allocMatrix refactor_count_alloc_matrix
#define R_alloc refactor_count_r_alloc
#define Rf_coerceVector refactor_count_coerce
#define INTEGER refactor_count_integer
#define REAL refactor_count_real
#define Rf_mkChar refactor_count_mkchar
#define Rf_setAttrib refactor_count_setattr
#define R_GetCCallable refactor_count_get_callable
namespace actual_candidate_counts {
#include SCEPTRE_COUNT_SOURCE
}
namespace actual_candidate_covariates {
#include SCEPTRE_COVARIATE_SOURCE
}
#undef Rf_allocVector
#undef Rf_allocMatrix
#undef R_alloc
#undef Rf_coerceVector
#undef INTEGER
#undef REAL
#undef Rf_mkChar
#undef Rf_setAttrib
#undef R_GetCCallable

// [[Rcpp::export]]
SEXP refactor_count_ownership_call(int kernel, Rcpp::List inputs, int boundary,
                                 int native, int error_kind, bool bounded_overflow = false) {
  observing = false;
  if (live_blocks != 0) Rcpp::stop("Previous probe retained native ownership.");
  requested_boundary = boundary; requested_native = native; native_error_kind = error_kind;
  refactor_count_overflow_limit = bounded_overflow ? 2 : std::numeric_limits<int>::max();
  boundary_count = native_count = max_live_blocks = failed_boundary = failed_native = unsafe_boundaries = 0;
  observing = true;
  Rcpp::RObject result;
  if (kernel == 0) {
    result = actual_candidate_counts::sparse_pair_counts(inputs[0], inputs[1], 4,
        inputs[2], inputs[3], inputs[4], inputs[5], false, inputs[6]);
  } else if (kernel == 1) {
    result = actual_candidate_counts::sparse_group_counts(inputs[0], inputs[1], 4,
        inputs[2], inputs[6]);
  } else if (kernel == 2) {
    result = actual_candidate_covariates::compute_cell_covariates_cpp(inputs[7], inputs[8],
        inputs[9], 4, 1, inputs[10], true, true);
  } else Rcpp::stop("Unknown count kernel.");
  // Track allocations through the function return and conversion to an R result too.
  Rcpp::RObject wrapped(Rcpp::wrap(result));
  observing = false;
  return wrapped;
}

// [[Rcpp::export]]
Rcpp::List refactor_count_ownership_state() {
  observing = false;
  return Rcpp::List::create(Rcpp::Named("boundary_count") = boundary_count,
      Rcpp::Named("native_count") = native_count, Rcpp::Named("live_blocks") = live_blocks,
      Rcpp::Named("max_live_blocks") = max_live_blocks, Rcpp::Named("failed_boundary") = failed_boundary,
      Rcpp::Named("failed_native") = failed_native, Rcpp::Named("unsafe_boundaries") = unsafe_boundaries,
      Rcpp::Named("boundary_kinds") = Rcpp::IntegerVector(observed_kinds, observed_kinds + boundary_count),
      Rcpp::Named("boundary_live_blocks") = Rcpp::IntegerVector(observed_live, observed_live + boundary_count));
}

// [[Rcpp::export]]
void refactor_count_ownership_negative_control() {
  observing = false;
  if (live_blocks != 0) Rcpp::stop("Previous probe retained native ownership.");
  boundary_count = native_count = max_live_blocks = failed_boundary = failed_native = unsafe_boundaries = 0;
  requested_boundary = requested_native = 0;
  observing = true;
  refactor_count_vector<int> owner(1);
  r_boundary(1); // Safely throws C++ instead of deliberately leaking at Rf_error.
}
