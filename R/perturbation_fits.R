# CRT perturbation models concern binary treatment membership, independently of
# the guide-assignment mixture method. Persistent state is O(design columns),
# never fitted probabilities, QR arrays, cell vectors or resampling draws.

#' Reconstruct a compact logistic perturbation fit's probabilities
#' @param fit Compact fitted_coefs and diagnostics from
#' fit_perturbation_model().
#' @param design Actual fitting design in its original row and column order.
#' @return Transient probability vector, using the original binomial inverse
#' link.
#' @noRd
reconstruct_perturbation_probabilities <- function(fit, design) {
    coefficients <- fit$fitted_coefs
    if (
        !is.matrix(design) ||
            !is.numeric(design) ||
            !is.numeric(coefficients) ||
            length(coefficients) != ncol(design) ||
            any(!is.finite(coefficients))
    ) {
        stop(
            "A compact perturbation fit requires finite coefficients aligned ",
            "to its design."
        )
    }
    if (
        !is.null(colnames(design)) &&
            !identical(names(coefficients), colnames(design))
    ) {
        stop(
            "Perturbation fit coefficients do not match the ordered design ",
            "columns."
        )
    }
    # Match glm.fit's arithmetic and C inverse-link, including its clipping.
    # plogis() and ad hoc clipping do not reproduce the fitted probabilities.
    eta <- as.vector(design %*% coefficients) + rep.int(0, nrow(design))
    stats::binomial()$linkinv(eta)
}

#' Fit treatment membership and propose a compact reusable logistic model
#' @param trt_idxs Treated fitting-row positions; duplicate positions mark once.
#' @param design Actual fitting design, including its ordered column names.
#' @return List with original fitted probabilities and compact fit or NULL.
#' @details Warnings are captured as diagnostics and continue to the caller once
#'   on an actual fit. Cache reuse reconstructs probabilities without replaying
#'   historical fitting warnings. Aliased/nonfinite coefficients and any failure
#'   of exact reconstruction decline caching without a second fit.
#' @noRd
fit_perturbation_model <- function(trt_idxs, design) {
    diagnostics <- new.env(parent = emptyenv())
    diagnostics$warnings <- character()
    logistic_fit <- withCallingHandlers(
        fit_perturbation_glm(trt_idxs, design),
        warning = function(condition) {
            diagnostics$warnings <- c(
                diagnostics$warnings,
                conditionMessage(condition)
            )
        }
    )
    probabilities <- logistic_fit$fitted.values
    coefficients <- logistic_fit$coefficients
    if (!is.numeric(coefficients) || any(!is.finite(coefficients))) {
        return(list(probabilities = probabilities, fit = NULL))
    }
    compact <- list(
        fitted_coefs = coefficients,
        diagnostics = list(
            rank = logistic_fit$rank,
            converged = logistic_fit$converged,
            boundary = logistic_fit$boundary,
            iterations = logistic_fit$iter,
            warnings = diagnostics$warnings
        )
    )
    reconstructed <- tryCatch(
        reconstruct_perturbation_probabilities(compact, design),
        error = function(condition) NULL
    )
    if (!identical(reconstructed, probabilities)) {
        compact <- NULL
    }
    list(probabilities = probabilities, fit = compact)
}
