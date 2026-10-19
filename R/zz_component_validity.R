# CellView is transient: no counts are copied into it. The existing numerical
# drivers use dollar access; all fields still have typed S4 ownership.
#' Internal analysis cell selection
#'
#' Small temporary mappings between source, analysis and fitting cell positions.
#' Dollar access exposes typed fields to the internal numerical drivers.
#' @name CellView-class
#' @aliases $,CellView-method
#' @keywords internal
NULL

setMethod("$", "CellView", function(x, name) slot(x, name))
setValidity("CellView", function(object) {
    if (length(object@source_cells) != length(object@fit_to_analysis) ||
        anyNA(c(object@source_cells, object@fit_to_analysis, object@trt_fit_pos)) ||
        any(object@source_cells < 1L) || any(object@fit_to_analysis < 1L) ||
        anyDuplicated(object@source_cells) || anyDuplicated(object@fit_to_analysis) ||
        any(object@trt_fit_pos < 1L | object@trt_fit_pos > length(object@source_cells)) ||
        !identical(object@n_trt, as.integer(length(object@trt_fit_pos)))) {
        return("CellView has inconsistent ordered source, analysis or treatment positions")
    }
    TRUE
})

setValidity("ScreenData", function(object) {
    for (store in list(object@response_matrix, object@grna_matrix)) {
        if (!is.null(store) && !methods::is(store, "dgRMatrix") && !methods::is(store, "BackedCountStore")) {
            return("Count stores must be normalized row-accessible Matrix or typed backed stores")
        }
    }
    if (is.null(object@response_matrix) || is.null(object@grna_matrix)) return(TRUE)
    response <- object@response_matrix
    guide <- object@grna_matrix
    if (is.null(dim(response)) || is.null(dim(guide))) return(TRUE) # closed backing
    if (ncol(response) != ncol(guide)) return("Modalities must have aligned cells")
    if (!is.null(colnames(response)) && !is.null(colnames(guide)) &&
        !identical(colnames(response), colnames(guide))) return("Modality cell order differs")
    if (nrow(object@covariate_data_frame) && nrow(object@covariate_data_frame) != ncol(response)) {
        return("Cell covariates must match the source cells")
    }
    TRUE
})

setValidity("sceptre_object", function(object) {
    if (!identical(object@schema_version, 3L)) return("Unsupported sceptre schema version")
    valid_spec <- validObject(object@spec, test = TRUE)
    if (!isTRUE(valid_spec)) return(valid_spec)
    valid_fit <- validObject(object@fits, test = TRUE)
    if (!isTRUE(valid_fit)) return(valid_fit)
    if (length(object@fits@response_fits)) {
        ids <- rownames(object@data@response_matrix)
        if (any(!names(object@fits@response_fits) %in% ids)) return("Cached response identifiers must exist in ScreenData")
        if (!identical(object@fits@response_identity, response_fit_identity(object, object@fits@response_identity$cells))) {
            return("Response cache identity must match actual data and design inputs")
        }
    }
    valid <- validObject(object@data, test = TRUE)
    if (!isTRUE(valid)) return(valid)
    if (length(object@data@response_matrix) &&
        !is.null(dim(object@data@response_matrix))) {
        n <- ncol(object@data@response_matrix)
        cells <- object@preparation@analysis_cells
        if (anyNA(cells) || anyDuplicated(cells) || any(cells < 1L | cells > n)) {
            return("Prepared cells must be unique source-cell positions")
        }
        design <- object@preparation@design_matrix
        if (length(design) && nrow(design) != n) return("Design must match source cells")
        assignments <- object@preparation@guide_assignments
        if (length(assignments) && any(vapply(assignments, function(x) {
            !is.numeric(x) || anyNA(x) || any(x != as.integer(x)) ||
                any(x < 1L | x > n) || anyDuplicated(x) > 0L
        }, logical(1)))) return("Assignments must use canonical source-cell coordinates")
    }
    TRUE
})


# This is a homogeneous cache for one actual fitting scope. Each entry contains
# only the coefficient vector and scalar dispersion, not per-cell fit arrays.
validate_response_fit_entries <- function(entries, identity) {
    if (!length(entries)) return(TRUE)
    ids <- names(entries)
    if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
        return("Response fit identifiers must be named and unique")
    }
    if (!all(c("data", "cells", "columns", "design", "method") %in% names(identity)) ||
        !length(identity$cells) || anyNA(identity$cells) || anyDuplicated(identity$cells) ||
        any(identity$cells < 1L) || !length(identity$columns)) {
        return("Response fit cache requires an ordered fitting identity")
    }
    for (fit in entries) {
        if (!is.list(fit) || !is.numeric(fit$fitted_coefs) ||
            any(!is.finite(fit$fitted_coefs)) ||
            !identical(names(fit$fitted_coefs), identity$columns) ||
            !is.numeric(fit$theta) || length(fit$theta) != 1L ||
            !is.finite(fit$theta) || fit$theta <= 0) {
            return("Response fit coefficients or dispersion do not match their design")
        }
    }
    TRUE
}

validate_perturbation_fit_entries <- function(entries, identity) {
    if (!length(entries)) return(TRUE)
    ids <- names(entries)
    if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
        return("Perturbation fit identifiers must be named and unique")
    }
    if (!all(c("cells", "columns", "design", "method", "runtime") %in% names(identity)) ||
        !is.integer(identity$cells) || anyNA(identity$cells) || anyDuplicated(identity$cells) ||
        any(identity$cells < 1L) || !length(identity$columns) ||
        !identical(identity$method, perturbation_fit_method_id()) ||
        !is.list(identity$runtime)) {
        return("Perturbation fit cache requires an ordered analysis identity")
    }
    for (fit in entries) {
        if (!is.list(fit) || !setequal(names(fit), c("fitted_coefs", "diagnostics", "identity")) ||
            anyDuplicated(names(fit))) return("Perturbation fit records require compact model and identity fields")
        key <- fit$identity
        diagnostics <- fit$diagnostics
        if (!is.list(diagnostics) ||
            !setequal(names(diagnostics), c("rank", "converged", "boundary", "iterations", "warnings")) ||
            anyDuplicated(names(diagnostics)) || !is.character(diagnostics$warnings) || anyNA(diagnostics$warnings) ||
            any(vapply(diagnostics[c("converged", "boundary")], function(x)
                !is.logical(x) || length(x) != 1L || anyNA(x), logical(1))) ||
            any(vapply(diagnostics[c("rank", "iterations")], function(x)
                !is.numeric(x) || length(x) != 1L || anyNA(x) || !is.finite(x) ||
                x < 0L || x > .Machine$integer.max || x != as.integer(x), logical(1)))) {
            return("Perturbation fit diagnostics are invalid")
        }
        if (!is.list(fit) || !is.numeric(fit$fitted_coefs) ||
            any(!is.finite(fit$fitted_coefs)) ||
            !identical(names(fit$fitted_coefs), identity$columns) ||
            !is.list(fit$diagnostics) ||
            !setequal(names(key), c("cells", "design", "treatment")) || anyDuplicated(names(key)) ||
            any(vapply(key, function(x) !is.character(x) || length(x) != 1L ||
                anyNA(x) || !grepl("^[[:xdigit:]]{64}$", x), logical(1)))) {
            return("Perturbation fit coefficients or compact input identity are invalid")
        }
    }
    TRUE
}

setValidity("FitCache", function(object) {
    response <- validate_response_fit_entries(object@response_fits, object@response_identity)
    if (!isTRUE(response)) return(response)
    validate_perturbation_fit_entries(object@perturbation_fits, object@perturbation_identity)
})
