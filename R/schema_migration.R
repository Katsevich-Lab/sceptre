# Serialized S4 attributes are inspected directly: the registered class may
# already have changed, so slotNames() cannot describe the physical old object.
#' @noRd
migrate_sceptre_object <- function(object) {
    if (!methods::is(object, "sceptre_object")) stop("Expected a `sceptre_object`.")
    saved <- attributes(object)
    version <- if (is.null(saved$schema_version)) 1L else saved$schema_version
    if (identical(version, 3L)) return(object)
    if (!version %in% c(1L, 2L)) stop("Unsupported sceptre saved-object schema version.")
    if (version == 1L && !all(c("response_matrix", "grna_matrix", "functs_called") %in% names(saved))) {
        stop("Unrecognized legacy sceptre saved-object layout.")
    }
    if (version == 2L && !all(c("data", "spec", "preparation", "fits", "results") %in% names(saved))) {
        stop("Unrecognized component sceptre saved-object layout.")
    }
    components <- if (version == 1L) {
        stats::setNames(rep(list(saved), 5L), c("data", "spec", "preparation", "fits", "results"))
    } else lapply(saved[c("data", "spec", "preparation", "fits", "results")], attributes)
    # Boundary aliases are intentionally confined here; current objects use
    # domain names and direct count stores rather than schema-one/two wrappers.
    legacy_preparation_names <- c(
        design_matrix = "covariate_matrix",
        nt_nonzero_matrix = "M_matrix",
        control_nonzero_totals = "n_nonzero_tot_vector",
        guide_assignments = "initial_grna_assignment_list",
        source_groups = "grna_assignments_raw",
        analysis_groups = "grna_assignments",
        guide_counts_per_cell = "grnas_per_cell",
        assignment_exclusions = "cells_w_zero_or_twoplus_grnas",
        analysis_cells = "cells_in_use")
    for (field in names(legacy_preparation_names)) {
        components$preparation[[field]] <- components$preparation[[legacy_preparation_names[[field]]]]
    }
    for (field in c("response_matrix", "grna_matrix")) {
        store <- components$data[[field]]
        if (is.list(store)) {
            if (length(store) > 1L) stop("Malformed legacy count-store wrapper.")
            store <- if (length(store)) store[[1L]] else NULL
            if (is.list(store) && !length(store)) store <- NULL
        }
        descriptor <- components$data$backend[[sub("_matrix$", "", field)]]
        components$data[field] <- list(as_count_store(store, descriptor))
    }
    out <- methods::new("sceptre_object")
    # Only current component payloads cross the boundary; removed counters,
    # execution flags and working budgets do not become a parallel live schema.
    for (owner in c("data", "spec", "preparation")) {
        component <- methods::slot(out, owner)
        fields <- setdiff(intersect(methods::slotNames(component), names(components[[owner]])), "stages")
        for (field in fields) methods::slot(component, field, check = FALSE) <- components[[owner]][[field]]
        methods::slot(out, owner) <- component
    }
    if (version == 1L) {
        out@data@grna_target_data_frame$grna_group <- NULL
    }
    # Prefer original schema-two requests when available. Earlier saved formats
    # only retain expanded left-join rows, whose ordered blocks can be reversed
    # without dropping repeated requests.
    request_certainty <- character()
    for (field in c("discovery_pairs", "positive_control_pairs")) {
        expanded <- components$spec[[field]]
        if (is.null(expanded)) expanded <- data.frame()
        methods::slot(out@preparation, field) <- expanded
        original <- components$spec$requested[[field]]
        recovered <- recover_legacy_requested_pairs(expanded, original,
            components$preparation$grna_target_data_frame,
            components$spec$grna_integration_strategy)
        methods::slot(out@spec, field) <- recovered$pairs
        request_certainty[field] <- recovered$certainty
    }
    metadata <- components$preparation$grna_target_data_frame
    groups <- metadata$grna_group
    if (is.null(groups) && length(out@spec@grna_integration_strategy)) {
        canonical <- out@data@grna_target_data_frame
        groups <- if (out@spec@grna_integration_strategy == "union") canonical$grna_target else
            ifelse(canonical$grna_target == "non-targeting", "non-targeting", canonical$grna_id)
    }
    if (!is.null(groups)) out@preparation@guide_groups <- as.character(groups)
    # Count row identifiers do not reveal the feature labels supplied at import
    # (for example, Ensembl IDs can have separate mitochondrial gene labels).
    # Schema one did not persist those labels; keep their provenance unknown.
    previous_results <- components$results
    old_records <- previous_results$records
    effective_fields <- c("B1", "B2", "B3", "calibration_group_size", "n_calibration_pairs", "negative_control_pairs")
    effective <- components$spec[intersect(effective_fields, names(components$spec))]
    analysis_fields <- c(run_calibration_check = "calibration_result", run_power_check = "power_result", run_discovery_analysis = "discovery_result")
    for (analysis in names(analysis_fields)) {
        table <- previous_results[[analysis_fields[[analysis]]]]
        old <- old_records[[analysis]]
        if (is.null(table)) table <- data.frame()
        if (nrow(table) || !is.null(old)) {
            out@results@records[[analysis]] <- migrate_historical_record(table, old, version, effective)
        }
    }
    fit_count <- length(components$fits$response_precomputations)
    out@results@import_metadata <- list(migrated_from_schema = as.integer(version),
        provenance = "unknown", discarded_unverified_fit_count = as.integer(fit_count),
        canonical_requests = request_certainty)
    # An import ID is a legacy reopen hint, never evidence of immutable content.
    legacy_id <- components$data$integer_id
    if (length(legacy_id) && is.null(count_store_descriptor(out@data@response_matrix))) {
        out@results@import_metadata$legacy_backend_id <- legacy_id
    }
    out <- initialize_screen_identity(out)
    if (!any(request_certainty == "unknown")) {
        out@preparation@stages <- migrate_preparation_stages(out, saved$functs_called)
    }
    validObject(out)
    out
}

# Recover multiplicity from the historical left-join traversal, rather than
# uniquing rows. Unrecognized layouts remain readable but require fresh setup.
#' @noRd
recover_legacy_requested_pairs <- function(expanded, original, metadata, strategy) {
    strip_group <- function(pairs) pairs[, setdiff(names(pairs), "grna_group"), drop = FALSE]
    if (is.data.frame(original)) {
        original <- strip_group(original)
        for (field in intersect(c("grna_target", "response_id"), names(original))) original[[field]] <- as.character(original[[field]])
        return(list(pairs = original, certainty = "saved_original"))
    }
    canonical <- strip_group(expanded)
    if (!nrow(expanded) || identical(strategy, "union") || !"grna_group" %in% names(expanded)) {
        return(list(pairs = canonical, certainty = "reconstructed"))
    }
    if (!strategy %in% c("singleton", "bonferroni") ||
        !all(c("grna_target", "grna_group") %in% names(metadata))) {
        return(list(pairs = canonical, certainty = "unknown"))
    }
    starts <- integer()
    position <- 1L
    while (position <= nrow(expanded)) {
        groups <- as.character(metadata$grna_group[metadata$grna_target == expanded$grna_target[position]])
        width <- length(groups)
        block <- seq.int(position, length.out = width)
        if (!width || max(block) > nrow(expanded) ||
            !identical(as.character(expanded$grna_group[block]), groups) ||
            !all(vapply(block, function(row) identical(as.list(canonical[row, , drop = FALSE]),
                as.list(canonical[position, , drop = FALSE])), logical(1)))) {
            return(list(pairs = canonical, certainty = "unknown"))
        }
        starts <- c(starts, position)
        position <- position + width
    }
    list(pairs = canonical[starts, , drop = FALSE], certainty = "reconstructed")
}

#' @noRd
migrate_historical_record <- function(table, old, version, effective = list()) {
    if (!is.list(old)) old <- list()
    history <- lapply(old$history, function(entry) {
        previous <- if (is.null(entry$result)) entry$table else entry$result
        if (is.null(previous)) previous <- data.frame()
        migrate_historical_record(previous, entry, version, effective)
    })
    list(table = table, identity = NULL,
        provenance = list(source = "legacy", schema_version = as.integer(version),
            certainty = "unknown", saved_provenance = old$provenance,
            effective_snapshot = effective),
        scope = if (is.null(old$scope)) "unknown" else old$scope,
        invalidated_by = "legacy provenance unknown", history = history)
}

# Flags are evidence only when their corresponding durable payload validates.
#' @noRd
migrate_preparation_stages <- function(object, flags) {
    prep <- object@preparation
    n <- nrow(object@data@covariate_data_frame)
    source_valid <- function(x) is.numeric(x) && !anyNA(x) &&
        !any(x != as.integer(x)) && !any(x < 1L | x > n) && !anyDuplicated(x)
    flag <- function(name) isTRUE(unname(flags[name]))
    object@preparation@stages <- list()
    add_stage <- function(stage) {
        store <- object@data@response_matrix
        if (!is.null(store) && (!methods::is(store, "BackedCountStore") || !is.null(store@handle))) {
            object <<- commit_preparation_stage(object, stage)
        } else {
            # Closed legacy stores cannot yet supply a live matrix to the
            # coordinator validator; payload checks below use saved cell axes.
            object@preparation@stages[[stage]] <<- list(identity = stage_input_identity(object, stage))
        }
    }
    design_ok <- flag("set_analysis_parameters") && n > 0L &&
        nrow(prep@design_matrix) == n && ncol(prep@design_matrix) > 0L &&
        length(prep@guide_groups) == nrow(object@data@grna_target_data_frame) &&
        !anyNA(prep@design_matrix)
    if (design_ok) add_stage("design")
    assignment_ok <- design_ok && flag("assign_grnas") &&
        length(prep@guide_assignments) > 0L &&
        all(c("grna_group_idxs", "indiv_nt_grna_idxs") %in% names(prep@source_groups)) &&
        all(vapply(prep@guide_assignments, source_valid, logical(1)))
    if (assignment_ok) add_stage("assignment")
    cells_ok <- assignment_ok && flag("run_qc") && length(prep@analysis_cells) > 0L &&
        source_valid(prep@analysis_cells) && is.list(prep@analysis_groups$grna_group_idxs)
    if (cells_ok) add_stage("cell_qc")
    pair_ok <- function(pairs, requested) {
        is.data.frame(pairs) && nrow(pairs) == nrow(requested) &&
            (!nrow(pairs) || (all(c("response_id", "grna_group", "n_nonzero_trt",
                "n_nonzero_cntrl", "pass_qc") %in% names(pairs)) &&
                is.logical(pairs$pass_qc) && !anyNA(pairs$pass_qc)))
    }
    if (cells_ok && pair_ok(prep@discovery_pairs_with_info, prep@discovery_pairs) &&
        pair_ok(prep@positive_control_pairs_with_info, prep@positive_control_pairs)) {
        add_stage("pair_qc")
    }
    object@preparation@stages
}

# Content identity is computed at import/migration, never once per response.
# Live backing handles are deliberately excluded: their identity is process-local.
#' @noRd
initialize_screen_identity <- function(object) {
    counts <- lapply(c("response_matrix", "grna_matrix"), function(field) {
        store <- methods::slot(object@data, field)
        descriptor <- count_store_descriptor(store)
        if (!is.null(descriptor)) return(list(content = descriptor$content_md5,
            dimensions = descriptor$dimensions, dimnames = descriptor$dimnames))
        store
    })
    object@data@identity <- digest::digest(list(counts = counts,
        covariates = object@data@covariate_data_frame,
        guides = object@data@grna_target_data_frame,
        low_moi = object@data@low_moi), algo = "sha256")
    object
}

#' @noRd
backed_store_descriptor <- function(store) {
    if (methods::is(store, "BackedCountStore")) return(store@descriptor)
    if (!inherits(store, "odm")) stop("Expected an open ondisc count matrix.")
    path <- normalizePath(attr(store, "h5_file"), mustWork = TRUE)
    list(path = path, content_md5 = unname(tools::md5sum(path)),
         dimensions = attr(store, "dimension"), dimnames = list(attr(store, "feature_ids"), NULL),
         integer_id = attr(store, "integer_id"))
}

#' @noRd
validate_backed_descriptor <- function(store, descriptor) {
    if (is.null(descriptor)) stop("Backing store has no anchored content descriptor; reimport the screen.")
    actual <- backed_store_descriptor(store)
    actual$content_md5 <- unname(tools::md5sum(actual$path))
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
    for (field in c("response_matrix", "grna_matrix")) {
        store <- methods::slot(object@data, field)
        if (methods::is(store, "BackedCountStore")) {
            validate_backed_descriptor(store, store@descriptor)
        }
    }
    invisible(NULL)
}
