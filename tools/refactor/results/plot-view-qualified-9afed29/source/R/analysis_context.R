# Source cells index store columns; analysis positions index analysis_to_source;
# fit positions index one ordered view. All R positions here are one based.

#' Validate an ordered source-cell map
#' @param source_cells Numeric integer-valued source positions, without missing
#'   values or repetitions; order is significant.
#' @param n_source Number of source columns.
#' @return Integer map in caller order. Empty input is allowed here; workflow QC
#'   separately rejects an empty retained-cell set. Inputs/RNG/files are unchanged.
#' @noRd
prepare_cell_map <- function(source_cells, n_source) {
    if (!is.numeric(source_cells) || anyNA(source_cells) ||
        any(source_cells != as.integer(source_cells)) ||
        any(source_cells < 1L | source_cells > n_source) ||
        anyDuplicated(source_cells)) {
        stop("Analysis source cells must be unique integer positions within the response store.")
    }
    as.integer(source_cells)
}

#' Normalize raw source membership to unique analysis positions
#' @param source_members Source-column memberships; order/repetitions are ignored.
#' @param analysis_to_source Valid ordered source-cell map.
#' @return Integer analysis positions in analysis order, including `integer()` for
#'   no members. This set-membership rule deliberately differs from NT-pool order.
#' @noRd
source_membership_to_analysis <- function(source_members, analysis_to_source) {
    which(analysis_to_source %in% source_members)
}

#' Encode analysis memberships in the existing saved NT-local schema
#' @param nt_groups_analysis Named list of nonempty NT groups after QC, with
#'   one-based analysis positions. Supported low-MOI NT groups are disjoint.
#' @return List `all_nt_idxs` (ordered first-occurrence analysis pool) and
#'   `indiv_nt_grna_idxs` (named group positions within that pool). An empty list
#'   yields an empty integer pool/list; calibration validates active group count.
#' @details Explicit matching preserves the saved partition/order without
#'   cumulative out-of-pool positions for overlaps. Does not mutate input or RNG.
#' @noRd
encode_nt_fit_assignments <- function(nt_groups_analysis) {
    nt_to_analysis <- as.integer(unique(unlist(nt_groups_analysis, use.names = FALSE)))
    groups_fit <- lapply(nt_groups_analysis, match, table = nt_to_analysis)
    list(indiv_nt_grna_idxs = groups_fit, all_nt_idxs = nt_to_analysis)
}

#' Resolve a compound calibration group in its supplied coordinate view
#' @param group_id One valid compound NT ID, with constituent names joined by `&`.
#' @param nt_groups Named NT membership list, normally in analysis coordinates.
#' @param low_moi Logical; low-MOI groups concatenate, high-MOI groups use the
#'   unique union in first-occurrence order.
#' @return Unnamed memberships retaining the declared order/duplicate rule;
#'   no coordinate conversion, mutation, sampling or storage access occurs.
#' @noRd
calibration_group_membership <- function(group_id, nt_groups, low_moi) {
    nt_ids <- strsplit(group_id, "&", fixed = TRUE)[[1L]]
    members <- stats::setNames(unlist(nt_groups[nt_ids]), NULL)
    if (!low_moi) members <- unique(members)
    members
}

#' Resolve request identities once against validated response/group labels
#' @param pairs Table with `response_id` and `grna_group` columns in request order.
#' @param response_ids,group_ids Unique canonical lookup labels.
#' @return Data frame with one-based `pair_row_id`, `response_idx`, `group_idx`.
#'   Rows and repeated requests are retained; public validation owns duplicate/ID
#'   rejection. Missing labels are not validated here; no input is mutated.
#' @noRd
resolve_pair_rows <- function(pairs, response_ids, group_ids) {
    data.frame(
        pair_row_id = seq_len(nrow(pairs)),
        response_idx = match(pairs$response_id, response_ids),
        group_idx = match(pairs$grna_group, group_ids)
    )
}

#' Normalize saved assignments into transient analysis/request coordinates
#' @param cells_in_use,n_source Ordered source cells and source-column count.
#' @param grna_assignments Saved QC-filtered assignments. Treatment groups index
#'   analysis cells; NT groups are NT-local only when NT controls are selected.
#' @param pairs,response_ids Valid requested pairs and source-row labels.
#' @param calibration_check,control_group_complement,low_moi Logical role flags.
#' @return Ordinary list containing source/NT maps, canonical analysis groups,
#'   response/group labels, request rows and row splits by source response/group.
#'   Splits retain request order; calibration group IDs use first occurrence.
#' @details Decodes the saved NT-local schema once, retaining the ordered pool.
#'   No fits, draw banks, persistent state, RNG or files are changed.
#' @noRd
prepare_index_context <- function(
    cells_in_use, n_source, grna_assignments, pairs,
    response_ids, calibration_check, control_group_complement, low_moi
) {
    analysis_to_source <- prepare_cell_map(cells_in_use, n_source)
    nt_to_analysis <- if (control_group_complement) integer() else grna_assignments$all_nt_idxs
    nt_groups_analysis <- grna_assignments$indiv_nt_grna_idxs
    if (!control_group_complement) {
        nt_groups_analysis <- lapply(nt_groups_analysis, function(nt_fit_pos) {
            nt_to_analysis[nt_fit_pos]
        })
    }
    group_ids <- if (calibration_check) {
        unique(as.character(pairs$grna_group))
    } else names(grna_assignments$grna_group_idxs)
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

#' Obtain one requested treatment membership in analysis coordinates
#' @param index Prepared index context.
#' @param group_id One validated treatment or compound calibration group ID.
#' @return One-based analysis memberships, with calibration concatenation/union
#'   semantics above. Does not construct a fit/source view or mutate the index.
#' @noRd
group_treatment_analysis <- function(index, group_id) {
    if (index$calibration_check) {
        calibration_group_membership(group_id, index$nt_groups_analysis, index$low_moi)
    } else {
        index$groups_analysis[[group_id]]
    }
}

#' Compose one common or guide-specific fit-cell view
#' @param index Prepared index context, with all positions one based.
#' @param group_id One group ID, or `NULL` for the common view without treatment.
#' @return List `fit_to_analysis`, `source_cells`, `trt_fit_pos`, `n_trt`.
#'   Complement uses all analysis cells; NT calibration
#'   uses the ordered NT pool; NT discovery uses treatment then NT. `NULL` under
#'   NT discovery gives all analysis cells so group workers can subset that view.
#' @details Order/repetitions follow memberships; controls occupy the remaining
#'   fit positions, without a separately materialized control map.
#'   For analysis->source `c(2,4,5,8,11)`, treatment `c(2,3)` and NT
#'   `c(4,1,5)`, NT discovery has source `c(4,5,8,2,11)` and treatment `1:2`.
#'   Views are transient, constructed at their sharing scope, never a bank of all
#'   guide-by-cell maps. No design, storage, fit, RNG or cache is changed.
#' @noRd
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
        n_trt = length(trt_fit_pos)
    )
}

#' Prepare a lean treatment lookup in the common fit-cell view
#' @param index Prepared index context, captured once for a run.
#' @return Closure taking a group ID and returning one-based `trt_idxs`, `n_trt`.
#'   NT calibration translates analysis memberships to ordered NT-pool positions;
#'   otherwise the common view uses analysis positions. NT discovery workers use
#'   their separate guide fit view, not this common-view lookup.
#' @details Resolves memberships only: no full source map/design per pair,
#'   all-guide view cache, RNG draw or mutation of the captured index.
#' @noRd
make_treatment_index_lookup <- function(index) {
    function(group_id) {
        trt_analysis_pos <- group_treatment_analysis(index, group_id)
        trt_fit_pos <- if (index$calibration_check && !index$control_group_complement) {
            match(trt_analysis_pos, index$nt_to_analysis)
        } else trt_analysis_pos
        list(trt_idxs = trt_fit_pos, n_trt = length(trt_fit_pos))
    }
}

#' Prepare one common matrix design and canonical analysis inputs from the facade
#' @param object Sceptre object with full source-row design and valid saved cache.
#' @param analysis One of `calibration_check`, `power_check`, `discovery_analysis`.
#' @param pairs Already validated eligible pair table in request order.
#' @return Transient list `response_store`, `index`, `common_view`, aligned matrix
#'   `design`, `pairs`, `analysis`, `shared_response_fit`, `test_options` and
#'   compact `response_fits`. Options contain staged budgets, their checked total,
#'   side and curve-fitting flag. NT discovery keeps guide-specific fits uncached.
#' @details Preserves matrix dimensions for accepted intercept-only formulas.
#'   The saved full source-row design/cache schema remains untouched; no full
#'   assay, fitted-vector cache, all-guide views, RNG or file writes are created.
#' @noRd
prepare_analysis_context <- function(object, analysis, pairs) {
    total_resamples <- checked_resampling_total(object@B1, object@B2, object@B3)
    store <- get_response_matrix(object)
    index <- prepare_index_context(object@cells_in_use, ncol(store), object@grna_assignments,
        pairs, rownames(store), identical(analysis, "calibration_check"),
        object@control_group_complement, object@low_moi)
    common_view <- make_fit_view(index)
    list(
        response_store = store,
        index = index,
        common_view = common_view,
        design = object@covariate_matrix[common_view$source_cells, , drop = FALSE],
        pairs = pairs,
        analysis = analysis,
        shared_response_fit = index$calibration_check || index$control_group_complement,
        test_options = list(B1 = object@B1, B2 = object@B2, B3 = object@B3,
            total_resamples = total_resamples, side_code = object@side_code,
            fit_parametric_curve = object@resampling_approximation == "skew_normal"),
        response_fits = object@response_precomputations
    )
}
