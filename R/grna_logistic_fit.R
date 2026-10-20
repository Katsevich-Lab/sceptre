# Internal helpers for the optional fast X | Z logistic nuisance fit.

prepare_grna_logistic_design <- function(covariate_matrix) {
    if (!is.matrix(covariate_matrix) || !is.numeric(covariate_matrix)) {
        stop("`covariate_matrix` must be a numeric matrix.", call. = FALSE)
    }
    if (nrow(covariate_matrix) < 2L || ncol(covariate_matrix) < 1L) {
        return(structure(
            list(
                n = nrow(covariate_matrix),
                p = ncol(covariate_matrix),
                valid = FALSE,
                reason = "unsupported_design_dimensions"
            ),
            class = c("sceptre_grna_logistic_design", "list")
        ))
    }
    if (!is.double(covariate_matrix)) {
        covariate_matrix <- matrix(
            as.double(covariate_matrix),
            nrow = nrow(covariate_matrix),
            ncol = ncol(covariate_matrix),
            dimnames = dimnames(covariate_matrix)
        )
    }
    prepared <- prepare_grna_logistic_design_cpp(
        Z = covariate_matrix,
        condition = TRUE,
        rank_tolerance = 1e-12
    )
    class(prepared) <- c("sceptre_grna_logistic_design", "list")
    prepared
}


validate_grna_logistic_prepared_design <- function(
    prepared_design,
    covariate_matrix
) {
    if (!inherits(prepared_design, "sceptre_grna_logistic_design")) {
        stop(
            "`prepared_design` must come from `prepare_grna_logistic_design()`.",
            call. = FALSE
        )
    }
    if (
        length(prepared_design$n) != 1L ||
            length(prepared_design$p) != 1L ||
            !identical(as.integer(prepared_design$n), nrow(covariate_matrix)) ||
            !identical(as.integer(prepared_design$p), ncol(covariate_matrix))
    ) {
        stop(
            "`prepared_design` is incompatible with `covariate_matrix` dimensions.",
            call. = FALSE
        )
    }
    invisible(NULL)
}


validate_grna_logistic_indices <- function(trt_idxs, n_cells) {
    if (
        !(is.integer(trt_idxs) || is.double(trt_idxs)) ||
            anyNA(trt_idxs) ||
            any(!is.finite(trt_idxs)) ||
            any(trt_idxs != trunc(trt_idxs))
    ) {
        stop(
            "`trt_idxs` must contain finite, integer-valued one-based indices.",
            call. = FALSE
        )
    }
    if (any(trt_idxs < 1) || any(trt_idxs > n_cells)) {
        stop("`trt_idxs` contains an out-of-range index.", call. = FALSE)
    }
    if (anyDuplicated(trt_idxs)) {
        stop("`trt_idxs` cannot contain duplicates.", call. = FALSE)
    }
    as.integer(trt_idxs)
}


make_grna_fit_details <- function(
    value,
    method_requested,
    method_used,
    fallback_reason,
    converged,
    iterations
) {
    list(
        value = value,
        method_requested = method_requested,
        method_used = method_used,
        fallback_reason = fallback_reason,
        converged = converged,
        iterations = iterations
    )
}
