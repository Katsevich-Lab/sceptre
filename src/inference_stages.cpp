#include "inference.h"
#include "sampling_safety.h"

namespace sceptre {
namespace inference {

// Use separate portions of the already generated draws: B1 from offset 0,
// B2 from B1, and B3 from B1+B2. Stage 1 uses add-one empirical p-values.
// If p <= .02, optionally try a skew-normal fit in stage 2; otherwise, or if
// that fit is declined, use the applicable empirical result. A reached fit
// requires at least 100 scores. Never generate new draws or pool stage samples.
// When B3 is zero, retain the previous stage's scores (B2 if a fit was attempted,
// otherwise B1), including any sorting/reversal by the fitter. The result keeps
// its own final score vector; full draws stay unchanged. Preserve missing
// skew-normal parameters and the existing tail/floor rules.
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
// Pass the result pointer and settings through R's error-protected conversion.
// This small stack record owns no C++ resources. The cleanup function below
// is responsible for deleting tail, including when R reports an error.
struct AssociationRecord {
  EffectEstimate effect;
  double z_orig;
  TailResult* tail;
  bool return_resampling_dist;
};

// Copy the finished result into R objects, allocating only the null scores from
// the stage actually returned. Do not reserve unused B2/B3 space for a stage-1
// result. Use the raw R API here: R errors can skip C++ destructors, so avoid
// Rcpp objects requiring cleanup, new C++ buffers and C++ exceptions. The caller
// has already installed R_UnwindProtect so an R error still frees the C++ result.
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

// Free the C++ result exactly once, after successful conversion or before R
// continues an error. Deleting this result must not allocate or throw an error.
void delete_association_tail(void* data, Rboolean) noexcept {
  AssociationRecord* record = static_cast<AssociationRecord*>(data);
  delete record->tail;
  record->tail = nullptr;
}
} // namespace

// Return the effect, score, p-value, stage, tail parameters and optional null scores.
// The caller protected an R continuation token before constructing C++ results.
// Use that token to install cleanup without an allocation gap that could leak
// tail if R errors. The conversion above copies only the selected final scores.
SEXP association_result_record(const EffectEstimate& effect, double z_orig,
                               TailResult* tail, bool return_resampling_dist,
                               SEXP continuation) {
  AssociationRecord record = {effect, z_orig, tail, return_resampling_dist};
  return R_UnwindProtect(association_record_leaf, &record,
                         delete_association_tail, &record, continuation);
}

} // namespace inference
} // namespace sceptre
