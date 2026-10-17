// Test cleanup of C++ results using controlled errors, without exhausting memory.
// Wrap relevant calls before including Rcpp templates, keeping the package's
// actual buffer lifetimes and calculation scopes unchanged.
// [[Rcpp::depends(BH)]]
#define Rf_allocVector probe_result_alloc_vector
#define Rf_allocMatrix probe_result_alloc_matrix
#define R_alloc probe_result_r_alloc
#define Rf_coerceVector probe_result_coerce
#define INTEGER probe_result_integer
#define REAL probe_result_real
#define Rf_mkChar probe_result_mkchar
#define Rf_setAttrib probe_result_setattr
#define R_GetCCallable probe_result_get_callable
#define R_MakeUnwindCont probe_result_make_cont
#define R_UnwindProtect probe_result_unwind
#define R_PreserveObject probe_result_preserve
#define Rf_protect probe_result_protect
#include <Rcpp.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <unordered_set>
#include <vector>
#include <boost/container_hash/hash.hpp>
#include <boost/math/distributions/skew_normal.hpp>
#undef Rf_allocVector
#undef Rf_allocMatrix
#undef R_alloc
#undef Rf_coerceVector
#undef INTEGER
#undef REAL
#undef Rf_mkChar
#undef Rf_setAttrib
#undef R_GetCCallable
#undef R_MakeUnwindCont
#undef R_UnwindProtect
#undef R_PreserveObject
#undef Rf_protect

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
SEXP R_MakeUnwindCont();
SEXP R_UnwindProtect(SEXP (*)(void*), void*, void (*)(void*, Rboolean), void*, SEXP);
void R_PreserveObject(SEXP);
SEXP Rf_protect(SEXP);
}

namespace {
bool observing = false;
int requested_boundary = 0, requested_native = 0;
int boundary_count = 0, native_count = 0, live_blocks = 0, live_payloads = 0;
int max_live_blocks = 0, max_live_payloads = 0, protected_depth = 0;
int failed_boundary = 0, failed_native = 0, unsafe_boundaries = 0, cleanup_calls = 0;
int kinds[4096], block_counts[4096], payload_counts[4096], protections[4096];
void* tracked[16384] = {};
std::size_t tracked_bytes[16384] = {}, live_bytes = 0, peak_bytes = 0;

void native_boundary() {
  if (!observing) return;
  if (protected_depth) throw std::runtime_error("Native allocation inside protected R leaf.");
  if (++native_count == requested_native) {
    failed_native = native_count;
    throw std::bad_alloc();
  }
}
void r_boundary(int kind) {
  if (!observing) return;
  if (boundary_count == 4096) throw std::runtime_error("Result probe capacity exceeded.");
  kinds[boundary_count] = kind;
  block_counts[boundary_count] = live_blocks;
  payload_counts[boundary_count] = live_payloads;
  protections[boundary_count++] = protected_depth;
  if ((live_blocks || live_payloads) && !protected_depth) {
    ++unsafe_boundaries;
    // Detect the deliberately unsafe case with a C++ exception, allowing its
    // destructors to free memory rather than leaking it through an R error.
    throw std::runtime_error("Unprotected R boundary with native result ownership.");
  }
  if (boundary_count == requested_boundary) {
    failed_boundary = boundary_count;
    requested_boundary = 0;
    Rf_error("Controlled result R-boundary failure.");
  }
}
SEXP precious_preserve(SEXP object) {
  r_boundary(9);
  typedef SEXP (*Preserve)(SEXP);
  static Preserve preserve = reinterpret_cast<Preserve>(R_GetCCallable("Rcpp", "Rcpp_precious_preserve"));
  return preserve(object);
}
struct ProtectedLeaf {
  SEXP (*fun)(void*);
  void* data;
  void (*cleanup)(void*, Rboolean);
  void* cleanup_data;
  bool owns_result;
};
SEXP enter_leaf(void* data) {
  ++protected_depth;
  ProtectedLeaf* leaf = static_cast<ProtectedLeaf*>(data);
  return leaf->fun(leaf->data);
}
void leave_leaf(void* data, Rboolean jump) noexcept {
  ProtectedLeaf* leaf = static_cast<ProtectedLeaf*>(data);
  if (leaf->owns_result) ++cleanup_calls;
  leaf->cleanup(leaf->cleanup_data, jump);
  --protected_depth;
}
}
extern "C" SEXP probe_result_alloc_vector(SEXPTYPE type, R_xlen_t length) {
  r_boundary(1); return Rf_allocVector(type, length);
}
extern "C" SEXP probe_result_alloc_matrix(SEXPTYPE type, int rows, int columns) {
  r_boundary(2); return Rf_allocMatrix(type, rows, columns);
}
extern "C" char* probe_result_r_alloc(R_SIZE_T n, int size) {
  r_boundary(3); return R_alloc(n, size);
}
extern "C" SEXP probe_result_coerce(SEXP object, SEXPTYPE type) {
  r_boundary(4); return Rf_coerceVector(object, type);
}
extern "C" int* probe_result_integer(SEXP object) { r_boundary(5); return INTEGER(object); }
extern "C" double* probe_result_real(SEXP object) { r_boundary(6); return REAL(object); }
extern "C" SEXP probe_result_mkchar(const char* value) { r_boundary(7); return Rf_mkChar(value); }
extern "C" SEXP probe_result_setattr(SEXP object, SEXP name, SEXP value) {
  r_boundary(8); return Rf_setAttrib(object, name, value);
}
extern "C" DL_FUNC probe_result_get_callable(const char* package, const char* name) {
  r_boundary(10);
  if (!std::strcmp(package, "Rcpp") && !std::strcmp(name, "Rcpp_precious_preserve"))
    return reinterpret_cast<DL_FUNC>(&precious_preserve);
  return R_GetCCallable(package, name);
}
extern "C" SEXP probe_result_make_cont() { r_boundary(11); return R_MakeUnwindCont(); }
extern "C" SEXP probe_result_unwind(SEXP (*fun)(void*), void* data,
    void (*cleanup)(void*, Rboolean), void* cleanup_data, SEXP token) {
  ProtectedLeaf leaf = {fun, data, cleanup, cleanup_data, live_blocks != 0 || live_payloads != 0};
  return R_UnwindProtect(enter_leaf, &leaf, leave_leaf, &leaf, token);
}
extern "C" void probe_result_preserve(SEXP value) { r_boundary(12); R_PreserveObject(value); }
extern "C" SEXP probe_result_protect(SEXP value) { r_boundary(13); return Rf_protect(value); }

template<class T> struct ResultAllocator {
  using value_type = T;
  ResultAllocator() noexcept = default;
  template<class U> ResultAllocator(const ResultAllocator<U>&) noexcept {}
  T* allocate(std::size_t count) {
    native_boundary();
    T* pointer = std::allocator<T>{}.allocate(count);
    if (observing) {
      int slot = 0;
      while (slot < 16384 && tracked[slot]) ++slot;
      if (slot == 16384) {
        std::allocator<T>{}.deallocate(pointer, count);
        throw std::runtime_error("Native probe capacity exceeded.");
      }
      tracked[slot] = pointer; tracked_bytes[slot] = count * sizeof(T);
      live_bytes += tracked_bytes[slot]; peak_bytes = std::max(peak_bytes, live_bytes);
      ++live_blocks; max_live_blocks = std::max(max_live_blocks, live_blocks);
    }
    return pointer;
  }
  void deallocate(T* pointer, std::size_t count) noexcept {
    for (int slot = 0; slot < 16384; ++slot) if (tracked[slot] == pointer) {
      tracked[slot] = nullptr; live_bytes -= tracked_bytes[slot]; --live_blocks; break;
    }
    std::allocator<T>{}.deallocate(pointer, count);
  }
  template<class U> bool operator==(const ResultAllocator<U>&) const noexcept { return true; }
  template<class U> bool operator!=(const ResultAllocator<U>&) const noexcept { return false; }
};
template<class T> using result_vector = std::vector<T, ResultAllocator<T>>;
template<class Key, class Hash> using result_set =
  std::unordered_set<Key, Hash, std::equal_to<Key>, ResultAllocator<Key>>;

// Add a test-only marker to the actual result structs to count live objects.
// Track copies and moves too, including temporary TailResult returns. Production
// code has no such marker or test option.
struct ResultPayloadLifetime {
  ResultPayloadLifetime() {
    native_boundary();
    if (observing) { ++live_payloads; max_live_payloads = std::max(max_live_payloads, live_payloads); }
  }
  ResultPayloadLifetime(const ResultPayloadLifetime&) : ResultPayloadLifetime() {}
  ResultPayloadLifetime(ResultPayloadLifetime&&) : ResultPayloadLifetime() {}
  ResultPayloadLifetime& operator=(const ResultPayloadLifetime&) { return *this; }
  ResultPayloadLifetime& operator=(ResultPayloadLifetime&&) { return *this; }
  ~ResultPayloadLifetime() { if (observing) --live_payloads; }
};

#define Rf_allocVector probe_result_alloc_vector
#define Rf_allocMatrix probe_result_alloc_matrix
#define R_alloc probe_result_r_alloc
#define Rf_coerceVector probe_result_coerce
#define INTEGER probe_result_integer
#define REAL probe_result_real
#define Rf_mkChar probe_result_mkchar
#define Rf_setAttrib probe_result_setattr
#define R_GetCCallable probe_result_get_callable
#define R_MakeUnwindCont probe_result_make_cont
#define R_UnwindProtect probe_result_unwind
#define R_PreserveObject probe_result_preserve
#define Rf_protect probe_result_protect
#include SCEPTRE_RESULT_INFERENCE
#include SCEPTRE_RESULT_SCORES
#include SCEPTRE_RESULT_TAIL
#include SCEPTRE_RESULT_STAGES
#include SCEPTRE_RESULT_CALIBRATION
#include SCEPTRE_RESULT_ENTRY
#undef Rf_allocVector
#undef Rf_allocMatrix
#undef R_alloc
#undef Rf_coerceVector
#undef INTEGER
#undef REAL
#undef Rf_mkChar
#undef Rf_setAttrib
#undef R_GetCCallable
#undef R_MakeUnwindCont
#undef R_UnwindProtect
#undef R_PreserveObject
#undef Rf_protect

// [[Rcpp::export]]
SEXP result_ownership_bank(Rcpp::List rows) {
  observing = false;
  Rcpp::XPtr<sceptre::inference::DrawBank> pointer(static_cast<sceptre::inference::DrawBank*>(nullptr), true);
  std::unique_ptr<sceptre::inference::DrawBank> bank(new sceptre::inference::DrawBank);
  bank->resize(rows.size());
  for (R_xlen_t row = 0; row < rows.size(); ++row) {
    Rcpp::IntegerVector source(rows[row]);
    (*bank)[row].assign(source.begin(), source.end());
  }
  R_SetExternalPtrAddr(pointer, bank.release());
  return pointer;
}

// [[Rcpp::export]]
SEXP result_ownership_call(int kernel, Rcpp::List input, SEXP bank, int boundary, int native) {
  observing = false;
  if (live_blocks || live_payloads || protected_depth) Rcpp::stop("Previous result ownership remains live.");
  requested_boundary = boundary; requested_native = native;
  boundary_count = native_count = max_live_blocks = max_live_payloads = failed_boundary =
    failed_native = unsafe_boundaries = cleanup_calls = 0;
  peak_bytes = live_bytes = 0;
  observing = true;
  Rcpp::RObject result;
  if (kernel == 0) {
    result = sample_calibration_groups(2, 2, 6, 4, 4, 3, .5);
  } else if (kernel == 1) {
    result = sample_calibration_pairs(input[0], input[1], input[2],
                                      input[3], input[4], input[5], input[6]);
  } else if (kernel == 2) {
    result = run_association_test(input["y"], input["mu"], input["a"], input["w"], input["D"],
        input["trt_idxs"], input["n_trt"], input["use_all_cells"], bank,
        input["B1"], input["B2"], input["B3"], input["fit_parametric_curve"],
        input["return_resampling_dist"], input["side_code"]);
  } else Rcpp::stop("Unknown result kernel.");
  Rcpp::RObject wrapped(Rcpp::wrap(result)); // final copy/preservation after cleanup
  observing = false;
  return wrapped;
}

// [[Rcpp::export]]
Rcpp::List result_ownership_state() {
  observing = false;
  return Rcpp::List::create(Rcpp::Named("boundary_count") = boundary_count,
      Rcpp::Named("native_count") = native_count, Rcpp::Named("live_blocks") = live_blocks,
      Rcpp::Named("live_payloads") = live_payloads, Rcpp::Named("max_live_blocks") = max_live_blocks,
      Rcpp::Named("max_live_payloads") = max_live_payloads, Rcpp::Named("cleanup_calls") = cleanup_calls,
      Rcpp::Named("protected_depth") = protected_depth, Rcpp::Named("failed_boundary") = failed_boundary,
      Rcpp::Named("failed_native") = failed_native, Rcpp::Named("unsafe_boundaries") = unsafe_boundaries,
      Rcpp::Named("peak_native_vector_bytes") = static_cast<double>(peak_bytes),
      Rcpp::Named("boundary_kinds") = Rcpp::IntegerVector(kinds, kinds + boundary_count),
      Rcpp::Named("boundary_blocks") = Rcpp::IntegerVector(block_counts, block_counts + boundary_count),
      Rcpp::Named("boundary_payloads") = Rcpp::IntegerVector(payload_counts, payload_counts + boundary_count),
      Rcpp::Named("boundary_protected") = Rcpp::IntegerVector(protections, protections + boundary_count));
}

// [[Rcpp::export]]
void result_ownership_negative_control() {
  observing = false;
  if (live_blocks || live_payloads) Rcpp::stop("Previous ownership remains live.");
  requested_boundary = requested_native = 0;
  boundary_count = native_count = max_live_blocks = max_live_payloads =
    failed_boundary = failed_native = unsafe_boundaries = cleanup_calls = 0;
  observing = true;
  ResultPayloadLifetime payload;
  result_vector<int> values(1);
  r_boundary(1); // C++ error unwinds safely; deliberately unsafe ordering detected.
}
