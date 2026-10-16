#' Build a finite, full-rank cell design without dropping rows
#'
#' @param covariate_data_frame Data frame with one row per source cell.
#' @param formula_object Formula passed to `model.matrix()`, including its
#'   intercept/contrast choices.
#' @return Numeric cell-by-design-column matrix in input row order. Errors on
#'   dropped rows, nonfinite entries or rank deficiency; does not modify inputs.
#' @noRd
convert_covariate_df_to_design_matrix <- function(
    covariate_data_frame,
    formula_object
) {
    global_cell_covariates_new <- stats::model.matrix(
        object = formula_object,
        data = covariate_data_frame
    )
    # verify that the global cell covariates are OK after transformation
    if (nrow(global_cell_covariates_new) != nrow(covariate_data_frame)) {
        stop(
            "Some rows of `covariate_data_frame` were lost after applying the `formula object`. This is likely due to NAs, which should be removed first."
        )
    }
    for (col_name in colnames(global_cell_covariates_new)) {
        vect <- global_cell_covariates_new[, col_name]
        if (any(is.infinite(vect)) || any(is.na(vect))) {
            stop(
                "The column `",
                col_name,
                "` of the `covariate_data_frame` after the `formula object` has been applied contains entries that are -Inf, Inf, or NA. Remove these entries."
            )
        }
    }
    # verify that matrix is not rank-deficient
    rank_matrix <- Matrix::rankMatrix(global_cell_covariates_new)
    if (rank_matrix != ncol(global_cell_covariates_new)) {
        stop(
            "The `formula_object` contains redundant information. This often occurs when one variable is `nested` inside another. For example, if the data frame contains a column `lane` with levels `lane_1`, `lane_2`, `lane_3`, and `lane_4` and a column `batch` with levels `batch_1` and `batch_2`, and if `lane_1` and `lane_2` are contained entirely within `batch_1` while `lane_3` and `lane_4` are contained entirely within `batch`2`, then `batch` is redundant. In this case `batch` should be removed from the formula object."
        )
    }
    return(global_cell_covariates_new)
}


#' Compute import-time cell summaries from in-memory count columns
#'
#' @param matrix_in Feature-by-cell R count matrix, converted to column access.
#' @param feature_names Character response labels aligned with rows; an empty
#'   vector disables mitochondrial summaries even when requested.
#' @param compute_p_mito,compute_max_feature Logical flags for optional summaries.
#' @return Data frame in column order with `n_nonzero` (values greater than 0.5)
#'   and `n_umis`, plus requested mito/max-feature columns. Max-feature ties take
#'   the first row; zero-total fractions retain native division results (NaN).
#' @noRd
compute_cell_covariates <- function(
    matrix_in,
    feature_names,
    compute_p_mito,
    compute_max_feature
) {
    # make response matrix column accessible
    matrix_in <- set_matrix_accessibility(
        matrix_in,
        make_row_accessible = FALSE
    )
    # get MT gene idxs
    if (compute_p_mito) {
        mt_gene_idxs <- grep(
            pattern = "^MT-",
            x = feature_names,
            ignore.case = TRUE
        )
        compute_p_mito <- length(mt_gene_idxs) >= 1
    } else {
        mt_gene_idxs <- integer()
    }
    # call the low-level function
    out <- storage_cell_covariates(
        matrix_in, mt_gene_idxs, compute_p_mito, compute_max_feature
    )
    ret <- data.frame(n_nonzero = out$n_nonzero, n_umis = out$n_umi)
    if (compute_p_mito) {
        ret$p_mito <- out$p_mito
    }
    if (compute_max_feature) {
        ret$frac_umis_max_feature <- out$frac_umis_max_feature
        ret$feature_w_max_expression <- rownames(matrix_in)[
            out$max_feature + 1L
        ]
    }
    return(ret)
}


#' Generate the native permutation bank for one prepared analysis
#'
#' @param grna_assignments Post-QC named memberships in the saved analysis/NT-local
#'   representation; empty targeting groups are ignored for NT-control size ranges.
#' @param B Integer total draw budget, normally `B1 + B2 + B3`.
#' @param calibration_check,control_group_complement Logical analysis/view flags.
#' @param calibration_group_size Number of NT guides pooled per calibration group.
#' @param n_cells Number of retained analysis cells, used for complement controls.
#' @return External pointer owning a native bank of zero-based indices. Its
#'   interpretation is specific to the selected Fisher-Yates/hybrid fit view;
#'   these are not source-cell IDs. Samplers use fixed native seeds without
#'   consuming R random draws; the pointer must not be persisted in saved objects.
#' @details Complement-control discovery draws in all retained cells with length
#'   equal to the largest treatment group. Complement calibration instead bounds
#'   that length by the sum of the largest pooled NT groups, capped at n_cells
#'   because overlapping memberships cannot enlarge the cell universe. NT-only
#'   calibration uses that pooled-size bound in the ordered NT fit universe.
#'   NT-control discovery uses the hybrid sampler with N NT controls and the
#'   smallest/largest nonempty treatment sizes m/M: for m <= k <= M, each size-k
#'   prefix is valid in its N+k fit universe. Fixed-size consumers read only their treatment
#'   prefix, not every stored bank entry. Generate the full B bank once here;
#'   stage selection must not redraw or reseed its segments.
#' @noRd
get_synthetic_permutation_idxs <- function(
    grna_assignments,
    B,
    calibration_check,
    control_group_complement,
    calibration_group_size,
    n_cells
) {
    check_native_scalar_count(B, "B")
    check_native_scalar_count(n_cells, "n_cells")
    if (calibration_check) {
        check_native_scalar_count(calibration_group_size, "calibration_group_size", 1)
        active_nt_count <- sum(lengths(grna_assignments$indiv_nt_grna_idxs) > 0L)
        if (calibration_group_size > active_nt_count)
            stop("`calibration_group_size` exceeds the active non-targeting gRNA count after QC (", active_nt_count, ").")
    }
    if (calibration_check && !control_group_complement) {
        indiv_nt_sizes <- vapply(
            grna_assignments$indiv_nt_grna_idxs,
            length,
            FUN.VALUE = integer(1)
        ) |>
            sort(decreasing = TRUE)
        M <- sum(as.double(indiv_nt_sizes[seq_len(calibration_group_size)]))
        n_control_cells <- length(grna_assignments$all_nt_idxs)
        out <- fisher_yates_sampler(n_tot = n_control_cells, M = M, B = B)
    } else if (calibration_check && control_group_complement) {
        indiv_nt_sizes <- vapply(
            grna_assignments$indiv_nt_grna_idxs,
            length,
            FUN.VALUE = integer(1)
        ) |>
            sort(decreasing = TRUE)
        # Overlapping high-MOI memberships can sum beyond the cell count; any
        # union used by calibration is bounded by the number of cells.
        M <- min(sum(as.double(indiv_nt_sizes[seq_len(calibration_group_size)])), n_cells)
        out <- fisher_yates_sampler(n_tot = n_cells, M = M, B = B)
    } else if (!calibration_check && !control_group_complement) {
        grna_group_sizes <- vapply(
            grna_assignments$grna_group_idxs,
            length,
            FUN.VALUE = integer(1)
        )
        grna_group_sizes <- grna_group_sizes[grna_group_sizes != 0L]
        range_grna_group_sizes <- range(grna_group_sizes)
        n_control_cells <- length(grna_assignments$all_nt_idxs)
        out <- hybrid_fisher_iwor_sampler(
            N = n_control_cells,
            m = range_grna_group_sizes[1],
            M = range_grna_group_sizes[2],
            B = B
        )
    } else if (!calibration_check && control_group_complement) {
        max_cells_per_grna_group <- vapply(
            grna_assignments$grna_group_idxs,
            length,
            FUN.VALUE = integer(1)
        ) |>
            max()
        out <- fisher_yates_sampler(
            n_tot = n_cells,
            M = max_cells_per_grna_group,
            B = B
        )
    }

    return(out)
}


#' Add integration-specific group IDs to metadata and requested pair tables
#'
#' @param sceptre_object Parameterized object with guide/target metadata and pair
#'   tables whose rows may repeat.
#' @return Updated object. Union uses target IDs without expanding rows; singleton
#'   and Bonferroni expand each pair through a many-to-many target/guide join.
#'   Guide order follows metadata, repeated pair requests are retained, and NT
#'   guides keep the reserved `non-targeting` group label.
#' @noRd
update_dfs_based_on_grouping_strategy <- function(sceptre_object) {
    grouping_strategy <- sceptre_object@grna_integration_strategy
    if (grouping_strategy == "union") {
        sceptre_object@grna_target_data_frame$grna_group <- sceptre_object@grna_target_data_frame$grna_target
        sceptre_object@discovery_pairs$grna_group <- sceptre_object@discovery_pairs$grna_target
        sceptre_object@positive_control_pairs$grna_group <- sceptre_object@positive_control_pairs$grna_target
    }
    if (grouping_strategy %in% c("singleton", "bonferroni")) {
        grna_target_data_frame <- sceptre_object@grna_target_data_frame
        discovery_pairs <- sceptre_object@discovery_pairs
        positive_control_pairs <- sceptre_object@positive_control_pairs

        # function to update targeting pairs
        update_targeting_pairs <- function(
            targeting_pairs,
            updated_grna_target_data_frame
        ) {
            targeting_pairs$grna_group <- NULL
            dplyr::left_join(
                targeting_pairs,
                updated_grna_target_data_frame |>
                    dplyr::select(grna_target, grna_group),
                by = "grna_target",
                relationship = "many-to-many"
            )
        }

        # update each data frame, appending a grna group column
        updated_grna_target_data_frame <- grna_target_data_frame |>
            dplyr::mutate(
                grna_group = ifelse(
                    grna_target == "non-targeting",
                    "non-targeting",
                    grna_id
                )
            )
        sceptre_object@grna_target_data_frame <- updated_grna_target_data_frame
        sceptre_object@discovery_pairs <- update_targeting_pairs(
            discovery_pairs,
            updated_grna_target_data_frame
        )
        sceptre_object@positive_control_pairs <- update_targeting_pairs(
            positive_control_pairs,
            updated_grna_target_data_frame
        )
    }
    return(sceptre_object)
}
