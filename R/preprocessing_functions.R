#' Apply the model formula while keeping one design row per original cell
#' @param covariate_data_frame Covariates with one row per original cell.
#' @param formula_object Formula for model.matrix(), including
#'   intercept/contrasts.
#' @return Numeric design matrix in input-cell order. Error if transformation
#'   drops
#'   cells, creates missing/infinite values or gives redundant design columns.
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


#' Calculate cell summaries from an in-memory feature-by-cell count matrix
#' @param matrix_in Count matrix, converted to column access for the
#'   calculation.
#' @param feature_names Labels matching feature rows. An empty vector disables
#'   mitochondrial summaries even when requested.
#' @param compute_p_mito,compute_max_feature Whether to include mitochondrial
#'   fractions and the largest feature's contribution/name.
#' @return Data frame in original cell order with n_nonzero (values greater than
#'   0.5), n_umis and the requested optional summaries. Largest-feature ties use
#'   the first row; a zero total can produce NaN fractions.
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
        matrix_in,
        mt_gene_idxs,
        compute_p_mito,
        compute_max_feature
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


#' Generate all permutation draws once for an analysis
#' @param grna_assignments Saved post-QC guide memberships. Targeting groups use
#'   positions among cells kept after QC; with NT controls, NT groups use
#'   positions within the NT-cell list. Ignore empty targeting groups when
#'   finding the range of treatment sizes.
#' @param B Total number of draws, normally B1 + B2 + B3.
#' @param calibration_check,control_group_complement Whether this is calibration
#'   and whether controls are all other kept cells.
#' @param calibration_group_size Number of NT guides combined per calibration
#'   group.
#' @param n_cells Number of cells kept after QC, used with complement controls.
#' @return An R external pointer to C++ memory containing the complete draws.
#'   Draws contain positions starting at 0 within the data used by a fit, not
#'   original matrix column numbers. The pointer is temporary and must not be
#'   saved in a sceptre object. C++ samplers use fixed seeds without advancing
#'   R's
#'   random-number stream.
#' @details Complement discovery samples among all kept cells, drawing enough
#'   positions for the largest treatment group. Calibration needs enough for a
#'   combination of NT groups: use the sum of the largest group sizes, capped at
#'   the available cell count when memberships overlap. NT calibration samples
#'   within the ordered NT-cell list.
#'
#'   NT discovery fits each guide's treatment cells followed by N NT cells. The
#'   hybrid sampler supports treatment sizes m through M in one set of draws:
#'   the first k entries of each draw are valid positions among N+k cells.
#'   Entries after those k may belong to a larger fit and must not be read for
#'   the smaller one. Generate all B draws now and let test stages select their
#'   portions later, without resampling or resetting seeds.
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
        check_native_scalar_count(
            calibration_group_size,
            "calibration_group_size",
            1
        )
        active_nt_count <- sum(
            lengths(grna_assignments$indiv_nt_grna_idxs) > 0L
        )
        if (calibration_group_size > active_nt_count) {
            stop(
                "`calibration_group_size` exceeds the active non-targeting ",
                "gRNA count after QC (",
                active_nt_count,
                ")."
            )
        }
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
        M <- min(
            sum(as.double(indiv_nt_sizes[seq_len(calibration_group_size)])),
            n_cells
        )
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


#' Add the guide-group names used by the chosen integration strategy
#' @param sceptre_object Object with guide/target metadata and requested pair
#'   tables.
#' @return Updated object. Union uses target names without expanding pairs.
#'   Singleton and Bonferroni join each requested target to all matching guides,
#'   so multiple guides or repeated metadata matches can create multiple rows.
#'   Keep metadata guide order and repeated requests. NT guides retain the
#'   reserved group label non-targeting.
#' @noRd
update_dfs_based_on_grouping_strategy <- function(sceptre_object) {
    metadata <- sceptre_object@data@grna_target_data_frame
    strategy <- sceptre_object@spec@grna_integration_strategy
    groups <- if (strategy == "union") {
        metadata$grna_target
    } else {
        ifelse(
            metadata$grna_target == "non-targeting",
            "non-targeting",
            metadata$grna_id
        )
    }
    sceptre_object@preparation@guide_groups <- groups
    compile <- function(requested) {
        requested$grna_group <- NULL
        if (strategy == "union") {
            requested$grna_group <- requested$grna_target
            return(requested)
        }
        mapping <- data.frame(
            grna_target = metadata$grna_target,
            grna_group = groups
        )
        dplyr::left_join(
            requested,
            mapping,
            by = "grna_target",
            relationship = "many-to-many"
        )
    }
    sceptre_object@preparation@discovery_pairs <- compile(
        sceptre_object@spec@discovery_pairs
    )
    sceptre_object@preparation@positive_control_pairs <- compile(
        sceptre_object@spec@positive_control_pairs
    )
    sceptre_object
}
