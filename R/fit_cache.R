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
    !identical(object@preparation@design_matrix, design) ||
        !identical(object@spec@control_group_complement, complement)
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
    object@data@low_moi &&
        (!identical(object@spec@grna_assignment_method, method) ||
            !identical(
                object@spec@grna_assignment_hyperparameters,
                hyperparameters
            ))
}

#' Invalidate the response and perturbation families according to their inputs
#' @param object Sceptre object; the caller must keep the returned object.
#' @param event Reason for checking: `parameters`, `assignments` or `analysis_cells`.
#' @param changed TRUE when design, control selection or retained cells changed,
#'   or when a low-MOI assignment change can affect response-fitting cells.
#' @return The updated object. Any assignment rerun clears perturbation fits;
#'   response fits retain their existing low-MOI/design/cell invalidation policy.
#'   At lookup, each family also verifies its actual fitting-input identity.
#' @noRd
invalidate_fit_caches <- function(object, event, changed) {
    stopifnot(event %in% c("parameters", "assignments", "analysis_cells"))
    if (changed || identical(event, "assignments")) {
        object@fits@perturbation_fits <- list()
        object@fits@perturbation_identity <- list()
    }
    if (changed) {
        object@fits@response_fits <- list()
        object@fits@response_identity <- list()
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
    combined <- c(cache, updates)
    if (
        length(combined) &&
            (is.null(names(combined)) || anyDuplicated(names(combined)))
    ) {
        stop("Response fit identifiers must be unique.")
    }
    combined
}

#' Save response fits after an analysis finishes successfully
#' @param object Sceptre object; the caller must keep the returned object.
#' @param cache Response fits valid for this object's current design and cells.
#' @return The object with `response_precomputations` replaced by `cache`.
#'   The public analysis functions call this after the tests finish successfully.
#'   It saves the supplied fits without fitting models or checking them again.
#' @noRd
commit_response_fits <- function(object, cache) {
    identity <- attr(cache, "fit_identity")
    attr(cache, "fit_identity") <- NULL
    candidate <- object@fits
    candidate@response_fits <- cache
    candidate@response_identity <- if (is.null(identity)) list() else identity
    validObject(candidate)
    if (
        length(cache) &&
            !identical(identity, response_fit_identity(object, identity$cells))
    ) {
        stop("Response fits do not match the actual fitting inputs.")
    }
    object@fits <- candidate
    validObject(object)
    object
}

# Hash only the actual design and ordered fitting cells; the immutable imported
# content identity is computed once, never once per response. Coefficients retain
# the current Poisson-fit/dispersion algorithm.
response_fit_identity <- function(object, source_cells) {
    list(
        data = object@data@identity,
        cells = as.integer(source_cells),
        columns = colnames(object@preparation@design_matrix),
        design = digest::digest(
            object@preparation@design_matrix[
                source_cells,
                ,
                drop = FALSE
            ],
            algo = "sha256"
        ),
        method = response_fit_method_id()
    )
}

# A content replacement invalidates every derived stage. Reopening a verified
# backing descriptor uses a separate path and does not call this mutator.
invalidate_screen_dependents <- function(object) {
    object <- invalidate_result_records(object, "screen content replaced")
    object@fits <- new("FitCache")
    object@preparation <- new("PreparationRecord")
    initialize_screen_identity(object)
}

# Cache maintenance invariants (PR #233, adapted to component ownership):
# response coefficients/dispersion share one actual ordered fitting scope;
# perturbation coefficients share the analysis context and each carry compact
# hashes of their exact ordered fitting cells/design/binary treatment indicator.
# Never repeat NT control vectors or save cell-length probabilities/CRT draws.
# Guide labels select candidates, not validity. Lookup must compare actual inputs.
# assign_grnas clears perturbation fits at either MOI; response invalidation keeps
# its existing low-MOI/design/cell policy. Calibration bypasses perturbation fits
# and preserves discovery entries. Commit only after all workers/results succeed.
# CRT draws remain freshly generated on every run. Fit warnings occur on actual
# fits; warm reuse retains diagnostics without re-emitting those fit warnings.
# Keep this note synchronized with invalidate_fit_caches(),
# commit_analysis_fits(), run_crt_analysis() and test-fit-cache-families.R.

perturbation_cache_identity <- function(object) {
    cells <- object@preparation@analysis_cells
    list(
        cells = cells,
        columns = colnames(object@preparation@design_matrix),
        design = digest::digest(
            object@preparation@design_matrix[cells, , drop = FALSE],
            algo = "sha256"
        ),
        method = perturbation_fit_method_id(),
        runtime = list(
            R = as.character(getRversion()),
            stats = as.character(utils::packageVersion("stats")),
            BLAS = unname(extSoftVersion()["BLAS"]),
            matprod = getOption("matprod")
        )
    )
}

perturbation_fit_identity <- function(source_cells, design, trt_idxs) {
    indicator <- integer(nrow(design))
    indicator[trt_idxs] <- 1L
    list(
        cells = digest::digest(as.integer(source_cells), algo = "sha256"),
        design = digest::digest(design, algo = "sha256"),
        treatment = digest::digest(indicator, algo = "sha256")
    )
}

lookup_perturbation_fit <- function(cache, group, identity) {
    fit <- cache[[group]]
    if (is.null(fit) || !identical(fit$identity, identity)) {
        return(NULL)
    }
    fit
}

merge_perturbation_fits <- function(cache, updates) {
    ids <- names(updates)
    if (
        length(updates) &&
            (is.null(ids) ||
                anyNA(ids) ||
                any(!nzchar(ids)) ||
                anyDuplicated(ids))
    ) {
        stop("Perturbation fit updates must be named and unique.")
    }
    for (id in ids) {
        if (is.null(updates[[id]])) {
            cache[id] <- NULL
        } else {
            cache[id] <- updates[id]
        }
    }
    cache
}

commit_perturbation_fits <- function(object, cache) {
    identity <- attr(cache, "fit_identity")
    attr(cache, "fit_identity") <- NULL
    if (
        length(cache) &&
            !identical(identity, perturbation_cache_identity(object))
    ) {
        stop(
            "Perturbation fits do not match the actual analysis design and cells."
        )
    }
    if (
        length(cache) &&
            !all(
                names(cache) %in%
                    names(object@preparation@analysis_groups$grna_group_idxs)
            )
    ) {
        stop(
            "Perturbation fit identifiers must belong to the prepared guide groups."
        )
    }
    object@fits@perturbation_fits <- cache
    object@fits@perturbation_identity <- if (is.null(identity)) {
        list()
    } else {
        identity
    }
    validObject(object)
    object
}

commit_analysis_fits <- function(object, output) {
    object <- commit_response_fits(object, output$response_precomputations)
    if (!is.null(output$perturbation_fits)) {
        object <- commit_perturbation_fits(object, output$perturbation_fits)
    }
    object
}
