# Response fits are compact coefficients/dispersion records in the existing S4
# list. View validity is owned by the supported setters, not by numerical fits.
# NT-control discovery does not cache its guide-specific fits; its stored fits
# belong to NT-only calibration. Complement analyses share the analysis view.
#' Detect supported parameter changes that invalidate response fits
#'
#' @param object Existing sceptre S4 facade.
#' @param design Proposed design matrix in analysis-cell order.
#' @param complement Proposed scalar logical control-group flag.
#' @return Scalar logical from strict `identical` comparisons, including design
#'   attributes and order. Does not mutate the object or inspect numerical fits.
#' @noRd
response_cache_parameters_changed <- function(object, design, complement) {
    !identical(object@covariate_matrix, design) ||
        !identical(object@control_group_complement, complement)
}

#' Detect assignment changes affecting the low-MOI response-fit view
#'
#' @param object Existing sceptre S4 facade.
#' @param method Proposed assignment-method string.
#' @param hyperparameters Proposed assignment-hyperparameter list.
#' @return Scalar logical. Low MOI compares method and hyperparameters strictly;
#'   high MOI returns false because its response-fit view is not assignment-filtered.
#'   Analysis-cell changes are invalidated separately. No assignment-fit cache exists.
#' @noRd
response_cache_assignments_changed <- function(object, method, hyperparameters) {
    object@low_moi &&
        (!identical(object@grna_assignment_method, method) ||
         !identical(object@grna_assignment_hyperparameters, hyperparameters))
}

#' Clear compact response fits after a supported facade change
#'
#' @param object Existing sceptre S4 facade; caller must retain the returned object.
#' @param event One of `parameters`, `assignments`, `analysis_cells`, checked here
#'   to keep the invalidation reason explicit at each caller.
#' @param changed Scalar logical already resolved by the supported setter/QC path.
#' @return Object with `response_precomputations` cleared iff `changed`; other
#'   slots are unchanged. No fit-view fingerprint or numerical check is inferred.
#' @noRd
invalidate_response_cache <- function(object, event, changed) {
    # Keep events explicit at callers, even where conservative invalidation
    # currently gives them the same implementation.
    stopifnot(event %in% c("parameters", "assignments", "analysis_cells"))
    if (changed) object@response_precomputations <- list()
    object
}

#' Look up a compact response fit in the existing saved-cache schema
#'
#' @param cache Named list of compact coefficient/dispersion records.
#' @param response_id One response ID used as the list key.
#' @return Matching record or `NULL` for a miss. Caller-owned facade invalidation
#'   guarantees the fit view is still valid; the response-ID key alone cannot do so.
#' @noRd
lookup_response_fit <- function(cache, response_id) cache[[response_id]]

#' Append worker cache misses in the existing insertion order
#'
#' @param cache Existing named list of compact response-fit records.
#' @param updates Named list of newly fitted responses, absent from `cache` and
#'   unique across workers; callers establish this invariant before fitting.
#' @return Concatenated list, preserving old records then update order. This is
#'   not an upsert: duplicate names are neither replaced nor deduplicated. No fitted
#'   vectors, designs, or per-response cell maps are persisted or inputs mutated.
#' @noRd
merge_response_fits <- function(cache, updates) {
    # Workers only return newly fitted responses. Preserve legacy insertion
    # order; no fitted vectors or per-response cell maps are persisted.
    c(cache, updates)
}

#' Commit completed analysis fits to the existing facade slot
#'
#' @param object Existing sceptre S4 facade; caller must retain the returned object.
#' @param cache Named list of compact response fits valid for the current fit view.
#' @return Object with only `response_precomputations` replaced. Analysis wrappers
#'   call this after successful execution; no validation, fitting or invalidation
#'   occurs here. No cache operation uses RNG or files.
#' @noRd
commit_response_fits <- function(object, cache) {
    object@response_precomputations <- cache
    object
}
