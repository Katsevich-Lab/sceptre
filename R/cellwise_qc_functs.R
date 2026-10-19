#' Select cells that pass cellwise QC and record the removal reasons
#'
#' Use the per-cell covariates, requested quantile ranges, mitochondrial cutoff,
#' explicit cell exclusions and precomputed zero/multiple-guide exclusions.
#' Keep values equal to either quantile endpoint; remove mitochondrial fractions
#' strictly above the cutoff.
#' @return Object with `cells_in_use` as original matrix columns in increasing
#'   order and `cell_removal_metrics` as counts by reason. Reasons can overlap,
#'   so their counts need not sum to the number of distinct removed cells.
#'   The caller checks for zero remaining cells, clears incompatible saved fits
#'   and updates guide assignments afterward. This step reads summaries rather
#'   than loading the whole count matrix.
#' @noRd
determine_cells_to_retain <- function(
    sceptre_object,
    response_n_umis_range,
    response_n_nonzero_range,
    p_mito_threshold,
    additional_cells_to_remove
) {
    # 1. remove additional cells specified by the user
    cells_to_exclude_user_specified <- additional_cells_to_remove
    n_cells_rm_user_specified <- length(cells_to_exclude_user_specified)

    # 2. compute the cells to retain based on n_umis and n_nonzero
    exclude_cells_by_clipping_counts <- function(v, percentile_range) {
        cutoffs <- stats::quantile(x = v, probs = percentile_range)
        which(v < cutoffs[1] | v > cutoffs[2])
    }
    cells_to_exclude_n_umis <- exclude_cells_by_clipping_counts(
        sceptre_object@data@covariate_data_frame$response_n_umis,
        response_n_umis_range
    )
    n_cells_rm_n_umis <- length(cells_to_exclude_n_umis)
    cells_to_exclude_n_responses <- exclude_cells_by_clipping_counts(
        sceptre_object@data@covariate_data_frame$response_n_nonzero,
        response_n_nonzero_range
    )
    n_cells_rm_n_responses <- length(cells_to_exclude_n_responses)

    # 3. compute cells to retain based on p_mito
    if ("response_p_mito" %in% colnames(sceptre_object@data@covariate_data_frame)) {
        cells_to_exclude_p_mito <- which(
            sceptre_object@data@covariate_data_frame$response_p_mito >
                p_mito_threshold
        )
        n_cells_rm_p_mito <- length(cells_to_exclude_p_mito)
    } else {
        cells_to_exclude_p_mito <- integer()
        n_cells_rm_p_mito <- 0L
    }

    # 4. compute cells to retain based on cells containing multiple grnas (if in
    # low MOI)
    cells_to_exclude_zero_twoplus_grnas <- sceptre_object@preparation@assignment_exclusions
    n_cells_rm_zero_twoplus_grnas <- length(
        sceptre_object@preparation@assignment_exclusions
    )

    # 5. finally, determine the set of cells to retain, and update the
    # sceptre_object
    cells_to_exclude <- c(
        cells_to_exclude_n_umis,
        cells_to_exclude_n_responses,
        cells_to_exclude_p_mito,
        cells_to_exclude_zero_twoplus_grnas,
        cells_to_exclude_user_specified
    ) |>
        unique()
    n_cells <- ncol(get_response_matrix(sceptre_object))
    cells_to_retain <- if (length(cells_to_exclude) == 0L) {
        seq(1L, n_cells)
    } else {
        seq(1L, n_cells)[-cells_to_exclude]
    }
    sceptre_object@preparation@analysis_cells <- cells_to_retain
    n_cells_rm_total <- n_cells - length(cells_to_retain)
    cell_removal_metrics <- c(
        n_cells_rm_n_umis = n_cells_rm_n_umis,
        n_cells_rm_n_responses = n_cells_rm_n_responses,
        n_cells_rm_p_mito = n_cells_rm_p_mito,
        n_cells_rm_zero_twoplus_grnas = n_cells_rm_zero_twoplus_grnas,
        n_cells_rm_user_specified = n_cells_rm_user_specified,
        n_cells_rm_total = n_cells_rm_total
    )
    sceptre_object@preparation@cell_removal_metrics <- cell_removal_metrics
    return(sceptre_object)
}


#' Remove failed-QC cells from guide assignments and renumber the remaining cells
#' @param sceptre_object Object with raw assignments as original matrix columns
#'   and `cells_in_use` listing the columns kept after QC.
#' @return Object with updated `grna_assignments`. Group members become unique
#'   positions within `cells_in_use`, in that list's order. Remove NT guides with
#'   no remaining cells. With NT controls, also store each NT group's positions
#'   within the combined NT-cell list, as required by the saved-object format.
#'   The caller must keep the returned object.
#' @noRd
update_grna_assignments_given_qc <- function(sceptre_object) {
    # Read the saved original-cell assignments and the cells retained by QC.
    cells_in_use <- sceptre_object@preparation@analysis_cells
    grna_assignments_raw <- sceptre_object@preparation@source_groups
    n_cells <- ncol(get_response_matrix(sceptre_object))

    analysis_to_source <- prepare_cell_map(cells_in_use, n_cells)

    # For every guide group, remove discarded cells and replace original
    # column numbers with positions in cells_in_use.
    grna_group_idxs_new <- lapply(
        grna_assignments_raw$grna_group_idxs,
        function(v) {
            source_membership_to_analysis(v, analysis_to_source)
        }
    ) |>
        stats::setNames(names(grna_assignments_raw$grna_group_idxs))
    nt_idxs_new <- lapply(grna_assignments_raw$indiv_nt_grna_idxs, function(v) {
        source_membership_to_analysis(v, analysis_to_source)
    }) |>
        stats::setNames(names(grna_assignments_raw$indiv_nt_grna_idxs))
    # remove those nt grnas with 0 cells (after QC)
    nt_idxs_new <- nt_idxs_new[
        vapply(nt_idxs_new, length, FUN.VALUE = integer(1)) != 0L
    ]
    grna_assignments <- list(
        grna_group_idxs = grna_group_idxs_new,
        indiv_nt_grna_idxs = nt_idxs_new
    )

    # With NT controls, save each NT group's positions within the combined
    # NT-cell list, rather than within the list of all cells kept after QC.
    if (!sceptre_object@spec@control_group_complement) {
        l <- encode_nt_fit_assignments(
            grna_assignments$indiv_nt_grna_idxs
        )
        grna_assignments$indiv_nt_grna_idxs <- l$indiv_nt_grna_idxs
        grna_assignments$all_nt_idxs <- l$all_nt_idxs
    }

    # 4. update the sceptre_object and return
    sceptre_object@preparation@analysis_groups <- grna_assignments
    return(sceptre_object)
}
