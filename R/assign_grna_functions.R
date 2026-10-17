###########################################
# TOP-LEVEL FUNCTION: ASSIGN GRNAS TO CELLS
###########################################
#' Assign guides to cells, then combine assignments when appropriate
#' @param sceptre_object Sceptre object with assignment method and settings chosen.
#' @param print_progress,parallel,n_processors,log_dir Mixture-model worker and
#'   progress settings.
#' @return Object with initial assignments as original matrix column numbers and,
#'   for maximum assignment, cells flagged for exclusion. Standard runs also
#'   combine assignments into guide groups. Nextflow runs process only
#'   `elements_to_analyze` and combine groups after workers' outputs are gathered.
#'   Warn if no guides were assigned. Mixture fitting can use workers and logs.
#' @noRd
assign_grnas_to_cells <- function(
    sceptre_object,
    print_progress,
    parallel,
    n_processors,
    log_dir
) {
    # extract pieces from sceptre_object
    grna_matrix <- get_grna_matrix(sceptre_object)
    grna_assignment_method <- sceptre_object@grna_assignment_method
    cell_covariate_data_frame <- sceptre_object@covariate_data_frame
    grna_assignment_hyperparameters <- sceptre_object@grna_assignment_hyperparameters

    # determine the gRNAs to analyze
    if (sceptre_object@nf_pipeline) {
        grnas_in_use <- sceptre_object@elements_to_analyze
    } else {
        grnas_in_use <- determine_grnas_in_use(sceptre_object)
    }

    # assign grnas via the selected strategy; obtain the grna
    # assignments and cells containing multiple grnas
    if (grna_assignment_method == "mixture") {
        sceptre_object@initial_grna_assignment_list <- assign_grnas_to_cells_mixture(
            grna_matrix = grna_matrix,
            cell_covariate_data_frame = cell_covariate_data_frame,
            grna_assignment_hyperparameters = grna_assignment_hyperparameters,
            print_progress = print_progress,
            parallel = parallel,
            n_processors = n_processors,
            log_dir = log_dir,
            grna_ids = grnas_in_use
        )
    }
    if (grna_assignment_method == "thresholding") {
        sceptre_object@initial_grna_assignment_list <- assign_grnas_to_cells_thresholding(
            grna_matrix = grna_matrix,
            grna_assign_threshold = grna_assignment_hyperparameters$threshold,
            grna_ids = grnas_in_use
        )
    }
    if (grna_assignment_method == "maximum") {
        max_result <- assign_grnas_to_cells_maximum(
            max_grna = sceptre_object@import_grna_assignment_info$max_grna,
            max_grna_frac_umis = sceptre_object@import_grna_assignment_info$max_grna_frac_umis,
            umi_fraction_threshold = grna_assignment_hyperparameters$umi_fraction_threshold,
            grna_n_umis = cell_covariate_data_frame$grna_n_umis,
            min_grna_n_umis_threshold = grna_assignment_hyperparameters$min_grna_n_umis_threshold,
            grna_ids = grnas_in_use
        )
        sceptre_object@initial_grna_assignment_list <- max_result$initial_assignment_list
        # set cells w/ multiple gRNAs for max method
        sceptre_object@cells_w_zero_or_twoplus_grnas <- max_result$cells_w_zero_or_twoplus_grnas
    }
    # check that at least some gRNAs were assigned
    if (
        all(
            vapply(
                sceptre_object@initial_grna_assignment_list,
                length,
                FUN.VALUE = integer(1)
            ) ==
                0L
        )
    ) {
        warning(
            "No gRNA was assigned to any cell. Consider setting `method` to 'thresholding' and `threshold` to a small, positive number."
        )
    }

    # process the initial assignment list
    if (!sceptre_object@nf_pipeline) {
        sceptre_object <- process_initial_assignment_list(sceptre_object)
    }
    return(sceptre_object)
}


#######################
# THRESHOLDING FUNCTION
#######################
#' Assign a guide to cells whose counts reach the threshold
#' @param grna_matrix Guide-by-cell matrix prepared for row access.
#' @param grna_assign_threshold Count cutoff.
#' @param grna_ids Guide names, in desired output order.
#' @return Named list of original cell-column numbers, starting at 1, in storage
#'   order. `storage_threshold_row()` handles in-memory versus ondisc counts,
#'   including their different numeric/integer threshold conversion.
#' @noRd
assign_grnas_to_cells_thresholding <- function(
    grna_matrix,
    grna_assign_threshold,
    grna_ids
) {
    grna_matrix <- set_matrix_accessibility(
        grna_matrix,
        make_row_accessible = TRUE
    )
    initial_assignment_list <- lapply(grna_ids, function(grna_id) {
        storage_threshold_row(
            grna_matrix,
            storage_response_index(grna_matrix, grna_id),
            grna_assign_threshold
        )
    }) |>
        stats::setNames(grna_ids)

    return(initial_assignment_list)
}

#########
# MAXIMUM
#########
#' Assign each cell to its highest-count guide and flag ambiguous cells
#' @param max_grna Highest-count guide name for each original cell.
#' @param max_grna_frac_umis Fraction of that cell's guide UMIs from its top guide.
#' @param umi_fraction_threshold Flag cells whose top-guide fraction is at or below this.
#' @param grna_n_umis Total guide UMIs per original cell.
#' @param min_grna_n_umis_threshold Flag cells whose total is below this.
#' @param grna_ids Guide names, in desired output order.
#' @return Named cell assignments and the combined set of ambiguous/low-total
#'   cells, as original matrix columns starting at 1. Flagged cells still appear
#'   in the initial assignments; subsequent QC removes them.
#' @noRd
assign_grnas_to_cells_maximum <- function(
    max_grna,
    max_grna_frac_umis,
    umi_fraction_threshold,
    grna_n_umis,
    min_grna_n_umis_threshold,
    grna_ids
) {
    initial_assignment_list <- lapply(grna_ids, function(grna_id) {
        which(max_grna == grna_id)
    }) |>
        stats::setNames(grna_ids)
    cells_w_multiple_grnas <- which(
        max_grna_frac_umis <= umi_fraction_threshold
    )
    cells_w_zero_grnas <- which(grna_n_umis < min_grna_n_umis_threshold)
    cells_w_zero_or_twoplus_grnas <- union(
        cells_w_multiple_grnas,
        cells_w_zero_grnas
    )
    return(list(
        initial_assignment_list = initial_assignment_list,
        cells_w_zero_or_twoplus_grnas = cells_w_zero_or_twoplus_grnas
    ))
}


######################
# HELPER FUNCTIONS
######################
#' Combine individual guide assignments into the groups saved before QC
#' @param sceptre_object Object with gathered guide assignments and guide metadata.
#' @return Object with `grna_assignments_raw` and `mean_cells_per_grna`. Methods
#'   other than maximum also count guides per cell and flag low-MOI cells with
#'   zero or multiple guides. Target groups follow first-seen metadata order;
#'   their cells are combined without duplicates, in first-seen order. Individual
#'   NT guide lists keep their original names and cell order. All cell numbers
#'   still refer to original matrix columns; QC renumbers them later.
#' @details Nextflow calls this after collecting guide assignments from its tasks.
#'   The supplied membership vectors are left unchanged.
#' @noRd
process_initial_assignment_list <- function(sceptre_object) {
    # set variables
    initial_assignment_list <- sceptre_object@initial_grna_assignment_list
    grna_target_data_frame <- sceptre_object@grna_target_data_frame
    low_moi <- sceptre_object@low_moi
    n_cells <- ncol(get_grna_matrix(sceptre_object))
    maximum_assignment <- sceptre_object@grna_assignment_method == "maximum"
    # if not using maximum assignment, compute n grnas per cell and cells
    # with multiple grnas
    if (!maximum_assignment) {
        sceptre_object@grnas_per_cell <- compute_n_grnas_per_cell_vector(
            initial_assignment_list,
            n_cells
        )
        if (low_moi) {
            cells_w_multiple_grnas <- which(sceptre_object@grnas_per_cell >= 2L)
            cells_w_zero_grnas <- seq(1, n_cells)[
                -sort(unique(unlist(initial_assignment_list)))
            ]
            sceptre_object@cells_w_zero_or_twoplus_grnas <- union(
                cells_w_multiple_grnas,
                cells_w_zero_grnas
            )
        }
    }
    # pool together targeting gRNAs via the or operation
    targeting_grna_group_data_table <- grna_target_data_frame |>
        dplyr::filter(grna_group != "non-targeting") |>
        data.table::as.data.table()
    targeting_grna_groups <- targeting_grna_group_data_table$grna_group |>
        unique()
    grna_group_idxs <- lapply(
        targeting_grna_groups,
        function(targeting_grna_group) {
            curr_grna_ids <- targeting_grna_group_data_table[
                targeting_grna_group_data_table$grna_group ==
                    targeting_grna_group,
            ]$grna_id
            initial_assignment_list[curr_grna_ids] |>
                unlist() |>
                unique()
        }
    ) |>
        stats::setNames(targeting_grna_groups)
    # obtain the individual non-targeting grna idxs
    nontargeting_grna_ids <- grna_target_data_frame |>
        dplyr::filter(grna_group == "non-targeting") |>
        dplyr::pull(grna_id)
    indiv_nt_grna_idxs <- initial_assignment_list[nontargeting_grna_ids]
    # construct the grna_group_idxs list
    grna_assignments_raw <- list(
        grna_group_idxs = grna_group_idxs,
        indiv_nt_grna_idxs = indiv_nt_grna_idxs
    )
    # initialize output
    sceptre_object@grna_assignments_raw <- grna_assignments_raw
    # save mean cells per grna
    sceptre_object@mean_cells_per_grna <- vapply(
        sceptre_object@initial_grna_assignment_list,
        length,
        FUN.VALUE = integer(1)
    ) |>
        mean()
    return(sceptre_object)
}


#' Choose which guides to process, retaining their metadata order
#' @param sceptre_object Object with guide metadata, requested targets and pipeline mode.
#' @param restricted_grnas FALSE selects every metadata row. TRUE selects guides
#'   for requested targets plus NT guides. For a nuclear trial run, instead sample
#'   at most 30 targeting metadata rows, take their targets, then add NT guides.
#'   Thirty rows need not mean 30 different targets.
#' @return Guide names in metadata order, keeping repeated metadata rows.
#'   Nuclear trial sampling uses seed 4 and restores R's previous random state.
#'   Nextflow's `output_grna_info` uses this restricted mode for trial runs.
#' @noRd
determine_grnas_in_use <- function(sceptre_object, restricted_grnas = FALSE) {
    grna_target_data_frame <- sceptre_object@grna_target_data_frame
    if (restricted_grnas) {
        if (!sceptre_object@nuclear) {
            all_grna_targets <- unique(c(
                sceptre_object@positive_control_pairs$grna_target,
                sceptre_object@discovery_pairs$grna_target,
                "non-targeting"
            ))
        } else {
            all_grna_targets <- c(
                withr::with_seed(4, {
                    grna_target_data_frame |>
                        dplyr::filter(grna_target != "non-targeting") |>
                        dplyr::sample_n(min(dplyr::n(), 30)) |>
                        dplyr::pull(grna_target)
                }),
                "non-targeting"
            )
        }
        grnas_in_use <- dplyr::filter(
            grna_target_data_frame,
            grna_target %in% all_grna_targets
        ) |>
            dplyr::pull(grna_id)
    } else {
        grnas_in_use <- grna_target_data_frame$grna_id
    }
    return(grnas_in_use)
}
