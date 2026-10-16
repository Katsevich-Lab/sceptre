# Response fits are compact coefficients/dispersion records in the existing S4
# list. View validity is owned by the supported setters, not by numerical fits.
# NT-control discovery does not cache its guide-specific fits; its stored fits
# belong to NT-only calibration. Complement analyses share the analysis view.
response_cache_parameters_changed <- function(object, design, complement) {
    !identical(object@covariate_matrix, design) ||
        !identical(object@control_group_complement, complement)
}

response_cache_assignments_changed <- function(object, method, hyperparameters) {
    object@low_moi &&
        (!identical(object@grna_assignment_method, method) ||
         !identical(object@grna_assignment_hyperparameters, hyperparameters))
}

invalidate_response_cache <- function(object, event, changed) {
    # Keep events explicit at callers, even where conservative invalidation
    # currently gives them the same implementation.
    stopifnot(event %in% c("parameters", "assignments", "analysis_cells"))
    if (changed) object@response_precomputations <- list()
    object
}

lookup_response_fit <- function(cache, response_id) cache[[response_id]]

merge_response_fits <- function(cache, updates) {
    # Workers only return newly fitted responses. Preserve legacy insertion
    # order; no fitted vectors or per-response cell maps are persisted.
    c(cache, updates)
}

commit_response_fits <- function(object, cache) {
    object@response_precomputations <- cache
    object
}
