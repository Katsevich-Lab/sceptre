# Turn test outputs into result tables, combine guides and adjust p-values.
# Table calculations take tables and settings directly. The workflow helpers
# near the end read those settings from the saved object and return a table.
#' Add each C++ test result to its requested gene-guide row
#'
#' The caller supplies one result per row, in the same order. Add p-values,
#' fold changes and their standard errors. Output level 2 also adds the stage,
#' observed score and fitted tail parameters; level 3 adds resampled statistics,
#' filling missing columns when tests returned different numbers of draws.
#' Guide integration and multiple-testing correction happen later.
#' @noRd
pair_results_from_native <- function(
    curr_df,
    curr_response_result,
    output_amount
) {
    curr_df$p_value <- vapply(
        X = curr_response_result,
        FUN = function(l) l$p,
        FUN.VALUE = numeric(1)
    )
    fold_change <- vapply(
        curr_response_result,
        FUN = function(l) l$fc,
        FUN.VALUE = numeric(1)
    )
    curr_df$fold_change <- fold_change
    curr_df$se_fold_change <- vapply(
        curr_response_result,
        FUN = function(l) l$se,
        FUN.VALUE = numeric(1)
    )
    curr_df$log_2_fold_change <- log(fold_change) / log(2)
    if (output_amount >= 2L) {
        curr_df$stage <- vapply(
            curr_response_result,
            FUN = function(l) l$stage,
            FUN.VALUE = integer(1)
        )
        curr_df$z_orig <- vapply(
            curr_response_result,
            FUN = function(l) l$z_orig,
            FUN.VALUE = numeric(1)
        )
        curr_df$xi <- vapply(
            curr_response_result,
            FUN = function(l) l$sn_params[1L],
            FUN.VALUE = numeric(1)
        )
        curr_df$omega <- vapply(
            curr_response_result,
            FUN = function(l) l$sn_params[2L],
            FUN.VALUE = numeric(1)
        )
        curr_df$alpha <- vapply(
            curr_response_result,
            FUN = function(l) l$sn_params[3L],
            FUN.VALUE = numeric(1)
        )
    }
    if (output_amount >= 3L) {
        to_append <- lapply(curr_response_result, FUN = function(l) {
            m <- data.table::as.data.table(matrix(l$resampling_dist, nrow = 1))
            colnames(m) <- paste0("z_null_", seq(1L, ncol(m)))
            m
        }) |>
            data.table::rbindlist(fill = TRUE)
        curr_df <- cbind(curr_df, to_append)
    }
    return(curr_df)
}

#' Combine guide-level results according to the chosen integration strategy
#'
#' For union, guide groups already represent targets: rename the group column.
#' For singleton or Bonferroni, join guide names to their targets, keeping all
#' metadata matches, then sort by p-value and response gene.
#'
#' For Bonferroni outside calibration, retain the smallest-p guide per gene-target
#' pair and multiply its p-value by the number of guides that passed QC, capped at
#' 1. If all guides failed QC, keep the largest treatment/control counts and leave
#' the effect and p-value missing. Calibration keeps its separate guide rows.
#' This guide-level adjustment precedes the later adjustment across tested pairs.
#' @noRd
integrate_guides <- function(
    result,
    grna_integration_strategy,
    grna_target_data_frame,
    is_calibration_check = FALSE
) {
    if (grna_integration_strategy == "union") {
        new_result <- result |> dplyr::rename("grna_target" = "grna_group")
    }
    if (grna_integration_strategy %in% c("singleton", "bonferroni")) {
        grna_target_data_frame <- grna_target_data_frame |>
            dplyr::select(grna_id, grna_target)
        new_result <- result |>
            dplyr::rename("grna_id" = "grna_group") |>
            dplyr::left_join(grna_target_data_frame, by = "grna_id") |>
            dplyr::relocate(response_id, grna_id, grna_target)
        if (
            grna_integration_strategy == "bonferroni" && !is_calibration_check
        ) {
            new_result <- new_result |>
                dplyr::group_by(response_id, grna_target) |>
                dplyr::group_modify(.f = function(tbl, key) {
                    if (!any(tbl$pass_qc)) {
                        n_nonzero_trt_out <- max(tbl$n_nonzero_trt)
                        n_nonzero_cntrl_out <- max(tbl$n_nonzero_cntrl)
                        log_2_fold_change_out <- p_value_out <- NA
                        pass_qc_out <- FALSE
                    } else {
                        min_p_idx <- which.min(tbl$p_value)
                        n_nonzero_trt_out <- tbl$n_nonzero_trt[min_p_idx]
                        n_nonzero_cntrl_out <- tbl$n_nonzero_cntrl[min_p_idx]
                        log_2_fold_change_out <- tbl$log_2_fold_change[
                            min_p_idx
                        ]
                        p_value_out <- min(
                            sum(tbl$pass_qc) * tbl$p_value[min_p_idx],
                            1
                        )
                        pass_qc_out <- TRUE
                    }
                    data.frame(
                        p_value = p_value_out,
                        n_nonzero_trt = n_nonzero_trt_out,
                        n_nonzero_cntrl = n_nonzero_cntrl_out,
                        log_2_fold_change = log_2_fold_change_out,
                        pass_qc = pass_qc_out
                    )
                }) |>
                dplyr::ungroup()
        }
        data.table::setorderv(
            new_result,
            c("p_value", "response_id"),
            na.last = TRUE
        )
    }
    return(new_result)
}

#' Put failed-QC pairs back into the result table
#'
#' Append the failed-QC requests to tested pairs, leaving unavailable test fields
#' missing. Sort by p-value and response gene, with missing p-values last, and
#' return a data.table. Keep pair labels and repeated requests. This sorts the
#' output; it does not restore original request order or adjust p-values.
#' @noRd
assemble_pair_results <- function(raw_result, requested_pairs) {
    result <- data.table::rbindlist(
        list(raw_result, requested_pairs |> dplyr::filter(!pass_qc)),
        fill = TRUE
    )
    data.table::setorderv(result, c("p_value", "response_id"), na.last = TRUE)
    result
}

#' Combine guides and mark significant results across the complete set of tests
#'
#' After guide integration, call `p.adjust()` with the chosen method and mark
#' adjusted p-values strictly below `alpha` as significant. Apply this to all
#' collected results, not independently to each Nextflow worker's partial output.
#' The scalar significance cutoff must come from the function argument: output
#' levels 2 and 3 also have an unrelated tail-shape column named `alpha`.
#' @noRd
finalize_association_result <- function(
    result,
    integration_strategy,
    target_map,
    correction_method,
    alpha,
    is_calibration_check = FALSE
) {
    result <- integrate_guides(
        result,
        integration_strategy,
        target_map,
        is_calibration_check
    )
    # Use the function's alpha cutoff, not the table's unrelated tail-shape
    # column also named alpha. Compute here, outside mutate's column lookup.
    significant <- stats::p.adjust(result$p_value, method = correction_method) <
        alpha
    dplyr::mutate(result, significant = !!significant)
}

# Nextflow scripts call these internal function names directly. Read settings
# from the object, then process the collected results. Standard analyses adjust
# p-values after collection; nuclear tasks keep their partial outputs raw.
#' Finish collected calibration results using the object's settings
#'
#' Called both by the package and directly by Nextflow scripts. Combine guides and
#' adjust across the complete set of results, keeping calibration guide rows
#' rather than using the Bonferroni row reduction for power/discovery. Return the
#' table; the caller saves it into the object if needed.
#' @noRd
process_calibration_result <- function(result, sceptre_object) {
    finalize_association_result(
        result,
        sceptre_object@spec@grna_integration_strategy,
        sceptre_object@preparation@grna_target_data_frame,
        sceptre_object@spec@multiple_testing_method,
        sceptre_object@spec@multiple_testing_alpha,
        TRUE
    )
}

#' Finish collected discovery results using the object's settings
#'
#' Called both by the package and directly by Nextflow scripts. Combine guides and
#' adjust p-values across the complete set of results. Return the table without
#' rerunning tests or saving it into the object.
#' @noRd
process_discovery_result <- function(result, sceptre_object) {
    finalize_association_result(
        result,
        sceptre_object@spec@grna_integration_strategy,
        sceptre_object@preparation@grna_target_data_frame,
        sceptre_object@spec@multiple_testing_method,
        sceptre_object@spec@multiple_testing_alpha
    )
}

#' Combine guide rows without marking significance across all tests
#'
#' Read guide metadata and the integration strategy from the object. The explicit
#' calibration flag determines whether calibration-specific grouping rules apply.
#' Return the grouped table. Power uses this step without the across-pair
#' significance calculation; nuclear pipeline tasks keep their raw-output rules.
#' @noRd
apply_grouping_to_result <- function(
    result,
    sceptre_object,
    is_calibration_check = FALSE
) {
    integrate_guides(
        result,
        sceptre_object@spec@grna_integration_strategy,
        sceptre_object@preparation@grna_target_data_frame,
        is_calibration_check
    )
}
