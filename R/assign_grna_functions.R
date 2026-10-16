###########################################
# TOP-LEVEL FUNCTION: ASSIGN GRNAS TO CELLS
###########################################
#' Dispatch guide assignment and preserve the pipeline pod boundary.
#' @param sceptre_object Facade with resolved assignment method/hyperparameters.
#' @param print_progress,parallel,n_processors,log_dir Existing mixture worker options.
#' @return Facade with initial source-cell memberships and maximum-method exclusion
#'   cells. Standard execution also pools/groups memberships; nf_pipeline uses
#'   elements_to_analyze and defers pooling until gathered assignments are processed.
#'   Warns if all assignments are empty. Reads storage; mixture may log and run GC.
#'   Assignment/source order and existing RNG behavior are retained.
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
            sceptre_object = sceptre_object,
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
#' Assign selected guide rows using the storage-specific threshold kernel.
#' @param grna_matrix Supported guide-by-source-cell store, made row-accessible.
#' @param grna_assign_threshold Existing scalar count cutoff.
#' @param grna_ids Requested guide IDs in output order.
#' @return Named list of one-based source-cell positions in stored-column order.
#'   Retains sparse double versus backed integer-coercion semantics at the storage
#'   boundary; does not pool guides, mutate the caller or draw randomness.
#' @noRd
assign_grnas_to_cells_thresholding <- function(
    grna_matrix,
    grna_assign_threshold,
    grna_ids
) {
    grna_matrix <- set_matrix_accessibility(grna_matrix, make_row_accessible = TRUE)
    initial_assignment_list <- lapply(grna_ids, function(grna_id) {
        storage_threshold_row(
            grna_matrix, storage_response_index(grna_matrix, grna_id),
            grna_assign_threshold
        )
    }) |>
        stats::setNames(grna_ids)

    return(initial_assignment_list)
}

#########
# MAXIMUM
#########
#' Resolve maximum-guide assignments from the existing import summaries.
#' @param sceptre_object Facade with source-cell max_grna and max_grna_frac_umis.
#' @param umi_fraction_threshold Cells at or below this fraction are flagged.
#' @param grna_n_umis Source-cell guide totals.
#' @param min_grna_n_umis_threshold Cells below this total are flagged.
#' @param grna_ids Guide IDs in output order.
#' @return Named one-based source-cell memberships and union of ambiguous/low-total
#'   positions. A flagged cell remains in the initial maximum membership; later QC
#'   applies exclusions. No matrix scan, input mutation, RNG or files.
#' @noRd
assign_grnas_to_cells_maximum <- function(
    sceptre_object,
    umi_fraction_threshold,
    grna_n_umis,
    min_grna_n_umis_threshold,
    grna_ids
) {
    initial_assignment_list <- lapply(grna_ids, function(grna_id) {
        which(sceptre_object@import_grna_assignment_info$max_grna == grna_id)
    }) |>
        stats::setNames(grna_ids)
    cells_w_multiple_grnas <- which(
        sceptre_object@import_grna_assignment_info$max_grna_frac_umis <=
            umi_fraction_threshold
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
# AUXILLIARY FUNCTIONS
######################
#' Pool source-cell assignments and prepare the saved raw-assignment schema.
#' @param sceptre_object Facade with gathered initial assignments and guide metadata.
#' @return Facade with grna_assignments_raw and mean_cells_per_grna; nonmaximum
#'   methods also count guides per source cell and low-MOI empty/multiple cells.
#'   Target groups retain first-seen metadata order and first-occurrence OR union;
#'   individual NT guides retain their initial named membership order. No QC view
#'   conversion occurs. This separate boundary is used after real pipeline pods.
#'   Does not mutate input memberships, sample, fit models or write files.
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


#' Select guide IDs in metadata order for standard or restricted/trial execution.
#' @param sceptre_object Facade containing metadata, pair targets and nuclear flag.
#' @param restricted_grnas Logical; FALSE retains all metadata rows. TRUE retains
#'   requested pair targets plus NT, or nuclear trial targets sampled from at most
#'   30 targeting metadata rows (not 30 distinct targets), plus NT.
#' @return Guide IDs in original metadata order, preserving metadata multiplicity.
#'   Nuclear sampling uses seed 4 and restores caller RNG. Retained because actual
#'   pipeline output_grna_info passes trial through this restricted-mode boundary.
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
