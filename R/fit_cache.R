# Save only response coefficients and dispersion, not fitted values for every
# cell. A saved fit is reusable only while its design and selected cells remain
# valid; the parameter setters and QC steps clear it when those inputs change.
# With NT controls, discovery fits different cells for each guide, so those fits
# are not saved here. Saved NT fits belong to calibration on the NT cells.
# With complement controls, analyses can share fits on all cells kept after QC.
#' Check whether changed analysis settings make saved response fits unusable
#' @param object Current sceptre object.
#' @param design Proposed design matrix, with one row per original cell.
#' @param complement Proposed setting for whether to use complement controls.
#' @return TRUE if the design or control setting differs from the saved value.
#'   Compare the entire matrix, including attributes and row order, even if only
#'   some of its cells are used in a fit. This function reports a change; it does
#'   not remove fits from the object.
#' @noRd
response_cache_parameters_changed <- function(object, design, complement) {
    !identical(object@covariate_matrix, design) ||
        !identical(object@control_group_complement, complement)
}

#' Check whether new guide assignments could change the response-fitting cells
#' @param object Current sceptre object.
#' @param method Proposed guide-assignment method.
#' @param hyperparameters Proposed guide-assignment settings.
#' @return For low MOI, TRUE if the method or settings changed. At high MOI,
#'   guide assignments do not select response-fitting cells, so return FALSE.
#'   Changes to the cells kept after QC are checked separately.
#' @noRd
response_cache_assignments_changed <- function(
    object,
    method,
    hyperparameters
) {
    object@low_moi &&
        (!identical(object@grna_assignment_method, method) ||
            !identical(object@grna_assignment_hyperparameters, hyperparameters))
}

#' Remove saved response fits when their inputs have changed
#' @param object Sceptre object; the caller must keep the returned object.
#' @param event Reason for checking: `parameters`, `assignments` or `analysis_cells`.
#' @param changed TRUE when the caller has detected a change requiring new fits.
#' @return The object with `response_precomputations` emptied if `changed` is TRUE.
#'   Other fields are unchanged. This function relies on the caller to decide
#'   whether fits are still usable; it does not compare their inputs itself.
#' @noRd
invalidate_response_cache <- function(object, event, changed) {
    # Name the reason at each call, even though all three reasons currently
    # clear the same saved-fit list.
    stopifnot(event %in% c("parameters", "assignments", "analysis_cells"))
    if (changed) {
        object@response_precomputations <- list()
    }
    object
}

#' Find a saved response fit by gene name
#' @param cache Named list of response coefficients and dispersion estimates.
#' @param response_id Gene name used as the list key.
#' @return The saved fit, or NULL if it is absent. Finding the gene name alone
#'   does not prove the fit is reusable: parameter setters and QC must already
#'   have cleared fits whose design or cell selection changed.
#' @noRd
lookup_response_fit <- function(cache, response_id) cache[[response_id]]

#' Append newly computed response fits to the saved-fit list
#' @param cache Existing named list of response fits.
#' @param updates New fits, named by gene. Callers must supply only previously
#'   absent genes and ensure that workers have not fitted the same gene twice.
#' @return One list containing the old fits followed by the new fits, in order.
#'   This is simple concatenation: it does not replace or remove duplicate names.
#'   Only coefficients and dispersion are stored, not cell-length fitted vectors.
#' @noRd
merge_response_fits <- function(cache, updates) {
    # Each worker returns only newly fitted genes. Keep the old fits first,
    # followed by the new coefficients and dispersion estimates.
    c(cache, updates)
}

#' Save response fits after an analysis finishes successfully
#' @param object Sceptre object; the caller must keep the returned object.
#' @param cache Response fits valid for this object's current design and cells.
#' @return The object with `response_precomputations` replaced by `cache`.
#'   The public analysis functions call this after the tests finish successfully.
#'   It saves the supplied fits without fitting models or checking them again.
#' @noRd
commit_response_fits <- function(object, cache) {
    object@response_precomputations <- cache
    object
}
