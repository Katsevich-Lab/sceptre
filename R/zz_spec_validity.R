# Correction policy is shared by public admission and component validity.
# Alpha=1 retains its existing valid interpretation and finite empirical budget;
# zero is unsupported because the empirical resampling-budget rule divides by it.
correction_settings_errors <- function(alpha, method) {
    problems <- character()
    if (!is.numeric(alpha) || is.complex(alpha) || !is.null(dim(alpha)) || length(alpha) != 1L || anyNA(alpha) ||
        !is.finite(alpha) || alpha <= 0 || alpha > 1) {
        problems["multiple_testing_alpha"] <-
            "`multiple_testing_alpha` must be one finite numeric value in (0, 1]."
    }
    supported <- is.character(method) && is.null(dim(method)) && length(method) == 1L && !anyNA(method)
    if (supported) {
        supported <- !inherits(tryCatch(match.arg(method, stats::p.adjust.methods),
            error = function(condition) condition), "error")
    }
    if (!supported) {
        problems["multiple_testing_method"] <- paste0(
            "`multiple_testing_method` must be one supported stats::p.adjust method: ",
            paste(stats::p.adjust.methods, collapse = ", "), ".")
    }
    problems
}

validate_correction_settings <- function(alpha, method) {
    problems <- correction_settings_errors(alpha, method)
    if (length(problems)) {
        condition <- structure(list(message = unname(problems[[1L]]), call = NULL,
            argument = names(problems)[[1L]]),
            class = c("sceptre_argument_error", "error", "condition"))
        stop(condition)
    }
    invisible(NULL)
}

analysis_spec_is_empty <- function(object) {
    # An untouched constructor prototype has neither a formula nor any payload.
    fields <- setdiff(slotNames(object), "formula_object")
    !is.call(object@formula_object) &&
        all(vapply(fields, function(field) length(slot(object, field)) == 0L, logical(1)))
}

validate_analysis_spec <- function(object) {
    if (!analysis_spec_is_empty(object)) {
        validate_correction_settings(object@multiple_testing_alpha,
            object@multiple_testing_method)
    }
    invisible(NULL)
}

setValidity("AnalysisSpec", function(object) {
    if (analysis_spec_is_empty(object)) return(TRUE)
    problems <- correction_settings_errors(object@multiple_testing_alpha,
        object@multiple_testing_method)
    if (length(problems)) unname(problems) else TRUE
})
