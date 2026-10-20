# Empty selector retains schema-three default certificates on old saved objects.
validate_perturbation_fit_method <- function(method) {
    if (
        !is.character(method) ||
            !is.null(dim(method)) ||
            length(method) > 1L ||
            anyNA(method) ||
            (length(method) && !method %in% c("glm.fit", "fast_logistic"))
    ) {
        stop("`method` must be one of 'glm.fit' or 'fast_logistic'.")
    }
    invisible(NULL)
}

#' Select the experimental CRT logistic fitter
#'
#' The reference default is `glm.fit`. `fast_logistic` reuses Ziang Niu's
#' unpenalized Newton solver, with conditioned design, convergence checks and
#' reference fallback. It targets the same likelihood but is not bit identical.
#' This selector concerns treatment membership, not guide assignment or SPA.
#'
#' @param sceptre_object An imported or prepared sceptre object.
#' @param method Exactly `"glm.fit"` or `"fast_logistic"`.
#' @return An updated object. Changing the method invalidates current results
#' and logistic caches; response fits and prepared QC remain reusable.
#' @details The experimental method requires a leading intercept and a finite,
#' full-rank design. Unsupported designs, nonconvergence, numerical rank loss,
#' degenerate treatment and suspected separation fall back to the original
#' `glm.fit` once, retaining its warnings. Compact fit diagnostics record the
#' requested/used method and fallback reason. Warm reuse retains diagnostics
#' without replaying warnings. Calibration uses the existing permutation-based
#' reference pathway; this selector affects fitted CRT treatment probabilities
#' for power and discovery. Validate calibration and power on each screen before
#' interpreting discoveries. No SPA method is enabled by this selector.
#' Package distributed execution rejects this experimental selector. Pipeline
#' use additionally requires an adapter selector guard; use the reference method
#' for pipeline workflows until that companion path is qualified.
#' @export
set_perturbation_fit_method <- function(sceptre_object, method = "glm.fit") {
    validate_perturbation_fit_method(method)
    if (!length(method)) {
        stop("`method` must be one scalar string.")
    }
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    if (!identical(get_perturbation_fit_method(sceptre_object), method)) {
        sceptre_object <- invalidate_result_records(
            sceptre_object,
            "CRT logistic fitter changed"
        )
        sceptre_object@fits@perturbation_fits <- list()
        sceptre_object@fits@perturbation_identity <- list()
    }
    sceptre_object@spec@perturbation_fit_method <- if (method == "glm.fit") {
        character()
    } else {
        method
    }
    validObject(sceptre_object)
    sceptre_object
}

#' @rdname set_perturbation_fit_method
#' @export
get_perturbation_fit_method <- function(sceptre_object) {
    method <- attr(sceptre_object@spec, "perturbation_fit_method", exact = TRUE)
    if (is.null(method) || !length(method)) {
        "glm.fit"
    } else {
        validate_perturbation_fit_method(method)
        method
    }
}

# The prepared design is private and constructed from the exact execution
# design. It is never persisted or accepted from a caller at the public API.
fit_fast_perturbation_glm <- function(
    trt_idxs,
    design,
    prepared_design = NULL
) {
    if (is.null(prepared_design)) {
        prepared_design <- prepare_grna_logistic_design(design)
    }
    # Existing indicator semantics treat duplicate positions once.
    trt_idxs <- validate_grna_logistic_indices(unique(trt_idxs), nrow(design))
    reason <- prepared_design$reason
    candidate <- NULL
    if (isTRUE(prepared_design$valid)) {
        candidate <- fit_grna_logistic_prepared_cpp(prepared_design, trt_idxs)
        reason <- candidate$reason
    }
    if (!is.null(candidate) && isTRUE(candidate$converged)) {
        result <- fast_logistic_result(candidate, design)
        if (!is.null(result)) {
            return(result)
        }
        reason <- "invalid_reconstructed_probabilities"
    }
    fit <- fit_perturbation_glm(trt_idxs, design)
    fit$candidate_details <- list(
        method_requested = "fast_logistic",
        method_used = "glm.fit",
        fallback_reason = reason
    )
    fit
}


# Use the canonical inverse link for compact warm/cold reconstruction.
fast_logistic_result <- function(candidate, design) {
    coefs <- as.numeric(candidate$coefficients)
    names(coefs) <- colnames(design)
    compact <- list(fitted_coefs = coefs)
    # Canonical link arithmetic makes cold and warm candidates bit identical
    # to each other, not to the reference IRLS fit.
    probability <- reconstruct_perturbation_probabilities(compact, design)
    if (
        all(is.finite(probability)) &&
            all(probability > 0 & probability < 1)
    ) {
        return(list(
            coefficients = coefs,
            fitted.values = probability,
            rank = candidate$rank,
            converged = TRUE,
            boundary = FALSE,
            iter = candidate$iterations,
            candidate_details = list(
                method_requested = "fast_logistic",
                method_used = "fast_logistic",
                fallback_reason = NA_character_
            )
        ))
    }
    NULL
}
