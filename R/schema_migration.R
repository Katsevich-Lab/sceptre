# Saved objects are explicitly upgraded from physical serialized attributes.
# The registered S4 definition may already be newer than an object's attributes.
#' @noRd
migrate_sceptre_object <- function(object) {
    if (!methods::is(object, "sceptre_object")) {
        stop("Expected a `sceptre_object`.")
    }
    saved <- attributes(object)
    if (all(c("schema_version", "data", "spec", "preparation", "fits", "results") %in% names(saved))) {
        if (!identical(saved$schema_version, 2L)) {
            stop("Unsupported sceptre saved-object schema version.")
        }
        return(object)
    }
    owners <- component_slot_owners()
    if (!all(c("response_matrix", "grna_matrix", "functs_called") %in% names(saved))) {
        stop("Unrecognized legacy sceptre saved-object layout.")
    }
    out <- methods::new("sceptre_object")
    for (field in intersect(names(owners), names(saved))) {
        component <- methods::slot(out, owners[[field]])
        methods::slot(component, field, check = FALSE) <- saved[[field]]
        methods::slot(out, owners[[field]]) <- component
    }
    for (field in intersect(c("last_function_called", "functs_called", "nuclear",
                             "nf_pipeline", "elements_to_analyze"), names(saved))) {
        methods::slot(out, field) <- saved[[field]]
    }
    out@preparation@grna_target_data_frame <- out@data@grna_target_data_frame
    # Flat legacy fields record the last effective settings, not necessarily
    # the user's original request (calibration can change the working budget).
    effective_names <- intersect(names(owners)[owners == "spec"], names(saved))
    out@spec@requested <- list(source = "legacy", certainty = "unknown",
        effective_snapshot = stats::setNames(lapply(effective_names,
            function(field) saved[[field]]), effective_names))
    out@results@records$legacy <- list(
        schema_version = 1L, provenance = "unknown", status = "historical",
        workflow = saved$functs_called, last_function_called = saved$last_function_called,
        unverified_fits = saved$response_precomputations,
        # Preserve fields unknown to this release without creating a second authority.
        extra_attributes = saved[setdiff(names(saved), c(names(owners),
            "last_function_called", "functs_called", "nuclear", "nf_pipeline",
            "elements_to_analyze", "class"))])
    for (analysis in c("calibration", "power", "discovery")) {
        table <- methods::slot(out@results, paste0(analysis, "_result"))
        if (nrow(table)) {
            record_name <- c(calibration = "run_calibration_check",
                power = "run_power_check", discovery = "run_discovery_analysis")[[analysis]]
            out@results@records[[record_name]] <- list(status = "historical",
                complete = NA, provenance = list(source = "legacy", certainty = "unknown"))
        }
    }
    out@fits@response_precomputations <- list()
    out@fits@identity <- list()
    final <- intersect(c("run_calibration_check", "run_power_check", "run_discovery_analysis"),
                       names(out@functs_called))
    out@functs_called[final] <- FALSE
    out@preparation@stages <- list(
        assignment = list(complete = isTRUE(unname(saved$functs_called["assign_grnas"])), provenance = "legacy"),
        qc = list(complete = isTRUE(unname(saved$functs_called["run_qc"])), provenance = "legacy"))
    initialize_screen_identity(out)
}

# Content identity is computed at import/migration, never once per response.
# Live backing handles are deliberately excluded: their identity is process-local.
#' @noRd
initialize_screen_identity <- function(object) {
    counts <- lapply(c("response_matrix", "grna_matrix"), function(field) {
        container <- methods::slot(object@data, field)
        if (!length(container)) return(NULL)
        store <- container[[1L]]
        if (methods::is(store, "odm")) {
            list(content = unname(tools::md5sum(store@h5_file)),
                 dimensions = dim(store), dimnames = dimnames(store))
        } else store
    })
    object@data@backend <- list()
    for (modality in c("response", "grna")) {
        container <- methods::slot(object@data, paste0(modality, "_matrix"))
        if (length(container) && methods::is(container[[1L]], "odm")) {
            object@data@backend[[modality]] <- backed_store_descriptor(container[[1L]])
        }
    }
    object@data@identity <- digest::digest(list(counts = counts,
        covariates = object@data@covariate_data_frame,
        guides = object@data@grna_target_data_frame,
        low_moi = object@data@low_moi), algo = "sha256")
    object
}

#' @noRd
backed_store_descriptor <- function(store) {
    if (!methods::is(store, "odm")) stop("Expected an open ondisc count matrix.")
    list(path = normalizePath(store@h5_file, mustWork = TRUE),
         content_md5 = unname(tools::md5sum(store@h5_file)),
         dimensions = dim(store), dimnames = dimnames(store),
         integer_id = store@integer_id)
}

#' @noRd
validate_backed_descriptor <- function(store, descriptor) {
    if (is.null(descriptor)) stop("Backing store has no anchored content descriptor; reimport the screen.")
    actual <- backed_store_descriptor(store)
    keys <- c("content_md5", "dimensions", "dimnames", "integer_id")
    if (!identical(actual[keys], descriptor[keys])) {
        stop("Backing store content does not match the saved descriptor; cached analyses cannot be reused.")
    }
    invisible(NULL)
}

# Saving verifies the anchored descriptor without rebaselining old provenance.
# During analysis, backing files are immutable by contract: getters/reporting
# and per-gene numerical loops do not scan files for unsupported external edits.
validate_backing_for_save <- function(object) {
    if (!length(object@data@identity)) return(invisible(NULL))
    for (modality in c("response", "grna")) {
        container <- methods::slot(object@data, paste0(modality, "_matrix"))
        if (length(container) && methods::is(container[[1L]], "odm")) {
            validate_backed_descriptor(container[[1L]], object@data@backend[[modality]])
        } else if (!length(container) && !is.null(object@data@backend[[modality]])) {
            descriptor <- object@data@backend[[modality]]
            actual <- unname(tools::md5sum(descriptor$path))
            if (is.na(actual) || !identical(actual, descriptor$content_md5)) {
                stop("Backing store content does not match the saved descriptor; cached analyses cannot be reused.")
            }
        }
    }
    invisible(NULL)
}
