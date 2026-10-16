#include "inference.h"
#include "sampling_safety.h"

namespace sceptre {
namespace inference {

// Consume disjoint segments of the already-generated full bank: B1 at offset
// 0, B2 at B1, B3 at B1+B2. Stage 1 uses add-one empirical tails; p<=.02 reaches
// optional skew-normal stage 2. A declined/failed fit uses empirical stage 3.
// B3=0 deliberately retains the previous selected segment (B2 if fitting was
// attempted, otherwise B1), including its fitter-produced sorting/reversal.
// Only reached B2 needs the fitter's >=100 observations. No stage redraws or
// concatenates earlier samples. Output owns its last segment; bank stays const,
// missing SN parameters stay NA and existing sidedness/floor conventions stay.
// No RNG or file effects; selected-stage validation errors propagate unchanged.
TailResult evaluate_staged_tail(double z_orig,
                                const Rcpp::NumericVector& a,
                                const Rcpp::NumericVector& w,
                                const Rcpp::NumericMatrix& D,
                                int n_trt, bool use_all_cells,
                                const DrawBank& synthetic_idxs,
                                int B1, int B2, int B3,
                                bool fit_parametric_curve, int side_code) {
  sceptre::sampling::sum(sceptre::sampling::sum(B1, B2, "B1 + B2"), B3,
                         "B1 + B2 + B3");
  const double P_THRESH = 0.02;
  bool sn_fit_used = false;
  TailResult out;
  out.sn_params = {{NA_REAL, NA_REAL, NA_REAL}};
  out.stage = 1;

  out.null_statistics = null_score_segment(a, w, D, 0, B1, n_trt,
                                          use_all_cells, synthetic_idxs);
  out.p = compute_empirical_p_value(out.null_statistics, z_orig, side_code);

  if (out.p <= P_THRESH) {
    if (fit_parametric_curve) {
      out.null_statistics = null_score_segment(a, w, D, B1, B2, n_trt,
                                              use_all_cells, synthetic_idxs);
      // The existing fitter may sort/reverse this owned vector. Preserve that
      // selected distribution, including when B3 is zero and the fit fails.
      std::vector<double> fit_sn_out =
        fit_and_evaluate_skew_normal(z_orig, out.null_statistics, side_code);
      out.p = fit_sn_out[3];
      sn_fit_used = out.p > -0.5;
      if (sn_fit_used) for (int i = 0; i < 3; i ++) out.sn_params[i] = fit_sn_out[i];
      out.stage = 2;
    }

    if (!fit_parametric_curve || !sn_fit_used) {
      if (B3 > 0) out.null_statistics = null_score_segment(a, w, D, B1 + B2,
                                                          B3, n_trt,
                                                          use_all_cells,
                                                          synthetic_idxs);
      out.p = compute_empirical_p_value(out.null_statistics, z_orig, side_code);
      out.stage = 3;
    }
  }
  return out;
}

namespace {
// A trivial stack context survives R's protected leaf; it owns no destructor.
// Only tail has native storage, deleted by the registered cleanup below.
struct AssociationRecord {
  EffectEstimate effect;
  double z_orig;
  TailResult* tail;
  bool return_resampling_dist;
};

// Exact reached-segment conversion, without Rcpp precious objects or native
// allocation/C++ throws. Do not allocate unused B2/B3 capacity for a stage1
// return. The caller registered cleanup before this leaf can allocate/error.
SEXP association_record_leaf(void* data) {
  const AssociationRecord* record = static_cast<AssociationRecord*>(data);
  const TailResult* tail = record->tail;
  const int fields = record->return_resampling_dist ? 7 : 6;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, fields));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, fields));
  const char* labels[] = {"p", "z_orig", "fc", "se", "stage", "sn_params", "resampling_dist"};
  const double scalars[] = {tail->p, record->z_orig, record->effect.fold_change,
                            record->effect.standard_error};
  for (int field = 0; field < fields; ++field)
    SET_STRING_ELT(names, field, Rf_mkChar(labels[field]));
  for (int field = 0; field < 4; ++field) {
    SET_VECTOR_ELT(out, field, Rf_allocVector(REALSXP, 1));
    REAL(VECTOR_ELT(out, field))[0] = scalars[field];
  }
  SET_VECTOR_ELT(out, 4, Rf_allocVector(INTSXP, 1));
  INTEGER(VECTOR_ELT(out, 4))[0] = tail->stage;
  SET_VECTOR_ELT(out, 5, Rf_allocVector(REALSXP, 3));
  double* sn_params = REAL(VECTOR_ELT(out, 5));
  for (int i = 0; i < 3; ++i) sn_params[i] = tail->sn_params[i];
  if (record->return_resampling_dist) {
    const R_xlen_t size = tail->null_statistics.size();
    SET_VECTOR_ELT(out, 6, Rf_allocVector(REALSXP, size));
    double* target = REAL(VECTOR_ELT(out, 6));
    for (R_xlen_t i = 0; i < size; ++i) target[i] = tail->null_statistics[i];
  }
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

// Runs exactly once after conversion or before R resumes its original jump.
// Destruction is allocation-free/noexcept, including an empty selected segment.
void delete_association_tail(void* data, Rboolean) noexcept {
  AssociationRecord* record = static_cast<AssociationRecord*>(data);
  delete record->tail;
  record->tail = nullptr;
}
} // namespace

// Preserve p/z/effect/SE/stage/SN and optional selected distribution, consuming
// native ownership at an exact-size R-only leaf. A precreated nonNULL protected
// token installs R's cleanup context without an allocating transfer window.
SEXP association_result_record(const EffectEstimate& effect, double z_orig,
                               TailResult* tail, bool return_resampling_dist,
                               SEXP continuation) {
  AssociationRecord record = {effect, z_orig, tail, return_resampling_dist};
  return R_UnwindProtect(association_record_leaf, &record,
                         delete_association_tail, &record, continuation);
}

} // namespace inference
} // namespace sceptre
