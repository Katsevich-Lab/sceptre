# Cell coordinates are R 1-based throughout this module. A source cell is a
# column in the store; an analysis position indexes analysis_to_source; a fit
# position indexes one ordered view. Contexts are transient ordinary lists.
prepare_cell_map <- function(source_cells, n_source) {
    if (!is.numeric(source_cells) || anyNA(source_cells) ||
        any(source_cells != as.integer(source_cells)) ||
        any(source_cells < 1L | source_cells > n_source) ||
        anyDuplicated(source_cells)) {
        stop("Analysis source cells must be unique integer positions within the response store.")
    }
    as.integer(source_cells)
}

source_membership_to_analysis <- function(source_members, analysis_to_source) {
    # The legacy mask operation returns unique memberships in analysis order,
    # even when raw assignments were reversed or contained duplicate cells.
    which(analysis_to_source %in% source_members)
}

legacy_nt_fit_assignments <- function(nt_groups_analysis) {
    nt_to_analysis <- unique(stats::setNames(unlist(nt_groups_analysis), NULL))
    sizes <- vapply(nt_groups_analysis, length, integer(1))
    stop <- cumsum(sizes)
    start <- c(0L, stop[-length(stop)]) + 1L
    groups_fit <- lapply(seq(1, length(nt_groups_analysis)), function(i) {
        seq(start[i], stop[i])
    }) |>
        stats::setNames(names(nt_groups_analysis))
    list(indiv_nt_grna_idxs = groups_fit, all_nt_idxs = nt_to_analysis)
}

calibration_group_membership <- function(group_id, nt_groups, low_moi) {
    nt_ids <- strsplit(group_id, "&", fixed = TRUE)[[1L]]
    members <- stats::setNames(unlist(nt_groups[nt_ids]), NULL)
    if (!low_moi) members <- unique(members)
    members
}

resolve_pair_rows <- function(pairs, response_ids, group_ids) {
    data.frame(
        pair_row_id = seq_len(nrow(pairs)),
        response_idx = match(pairs$response_id, response_ids),
        group_idx = match(pairs$grna_group, group_ids)
    )
}

prepare_index_context <- function(
    cells_in_use, n_source, grna_assignments, pairs,
    response_ids, calibration_check, control_group_complement, low_moi
) {
    analysis_to_source <- prepare_cell_map(cells_in_use, n_source)
    nt_to_analysis <- if (control_group_complement) integer() else grna_assignments$all_nt_idxs
    nt_groups_analysis <- grna_assignments$indiv_nt_grna_idxs
    if (!control_group_complement) {
        # The facade persists NT-local groups. Normalize once at the boundary;
        # retain the ordered NT pool instead of sorting the unique union.
        nt_groups_analysis <- lapply(nt_groups_analysis, function(nt_fit_pos) {
            nt_to_analysis[nt_fit_pos]
        })
    }
    group_ids <- if (calibration_check) unique(as.character(pairs$grna_group)) else names(grna_assignments$grna_group_idxs)
    pair_rows <- resolve_pair_rows(pairs, response_ids, group_ids)
    list(
        analysis_to_source = analysis_to_source,
        nt_to_analysis = nt_to_analysis,
        groups_analysis = grna_assignments$grna_group_idxs,
        nt_groups_analysis = nt_groups_analysis,
        pair_rows = pair_rows,
        response_ids = response_ids,
        group_ids = group_ids,
        rows_by_response = split(pair_rows$pair_row_id, pair_rows$response_idx),
        rows_by_group = split(pair_rows$pair_row_id, pair_rows$group_idx),
        calibration_check = calibration_check,
        control_group_complement = control_group_complement,
        low_moi = low_moi
    )
}

group_treatment_analysis <- function(index, group_id) {
    if (index$calibration_check) {
        calibration_group_membership(group_id, index$nt_groups_analysis, index$low_moi)
    } else {
        index$groups_analysis[[group_id]]
    }
}

make_fit_view <- function(index, group_id = NULL) {
    treatment <- if (is.null(group_id)) integer() else group_treatment_analysis(index, group_id)
    if (index$control_group_complement) {
        fit_to_analysis <- seq_along(index$analysis_to_source)
        trt_fit_pos <- treatment
    } else if (index$calibration_check) {
        fit_to_analysis <- index$nt_to_analysis
        trt_fit_pos <- match(treatment, fit_to_analysis)
    } else if (is.null(group_id)) {
        fit_to_analysis <- seq_along(index$analysis_to_source)
        trt_fit_pos <- integer()
    } else {
        fit_to_analysis <- c(treatment, index$nt_to_analysis)
        trt_fit_pos <- seq_along(treatment)
    }
    list(
        fit_to_analysis = fit_to_analysis,
        source_cells = index$analysis_to_source[fit_to_analysis],
        trt_fit_pos = trt_fit_pos,
        control_fit_pos = setdiff(seq_along(fit_to_analysis), trt_fit_pos),
        n_trt = length(trt_fit_pos)
    )
}

# Compatibility entry points still accept the historical assignment arguments.
# New drivers attach their prepared context to the one per-run index closure.
prepared_index_factory <- function(index) {
    function(group_id) {
        trt_analysis_pos <- group_treatment_analysis(index, group_id)
        trt_fit_pos <- if (index$calibration_check && !index$control_group_complement) {
            match(trt_analysis_pos, index$nt_to_analysis)
        } else trt_analysis_pos
        list(trt_idxs = trt_fit_pos, n_trt = length(trt_fit_pos))
    }
}

prepare_analysis_context <- function(object, analysis, pairs = NULL) {
    if (is.null(pairs)) {
        pairs <- switch(analysis,
            calibration_check = object@negative_control_pairs,
            power_check = object@positive_control_pairs_with_info,
            discovery_analysis = object@discovery_pairs_with_info,
            stop("Unknown analysis role: ", analysis)
        )
        if ("pass_qc" %in% names(pairs)) pairs <- pairs[pairs$pass_qc, ]
    }
    store <- get_response_matrix(object)
    list(
        response_store = store,
        index = prepare_index_context(object@cells_in_use, ncol(store), object@grna_assignments,
            pairs, rownames(store), identical(analysis, "calibration_check"),
            object@control_group_complement, object@low_moi),
        design = object@covariate_matrix,
        pairs = pairs,
        analysis = analysis,
        test_options = list(B1 = object@B1, B2 = object@B2, B3 = object@B3,
            side_code = object@side_code, resampling_approximation = object@resampling_approximation),
        caches = object@response_precomputations
    )
}

# The compatibility native inference entry point subtracts one internally for
# observed indices. Synthetic banks already contain zero-based indices and do
# not pass through this adapter.
native_observed_positions <- function(trt_fit_pos) as.integer(trt_fit_pos)

# Historical NT-union entry point retained for characterization/component use.
get_idx_vector_calibration_check <- function(
    curr_grna_group,
    indiv_nt_grna_idxs,
    low_moi
) {
    trt_idxs <- calibration_group_membership(curr_grna_group, indiv_nt_grna_idxs, low_moi)
    return(list(trt_idxs = trt_idxs, n_trt = length(trt_idxs)))
}

# Compatibility path for older direct worker callers. Current plans supply
# make_fit_view()'s composed source map, shared with their design matrix.
legacy_fit_source_cells <- function(analysis_to_source, nt_to_analysis, nt_only) {
    if (nt_only) analysis_to_source[nt_to_analysis] else analysis_to_source
}
