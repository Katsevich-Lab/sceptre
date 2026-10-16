# Ordinary raw-table operations; these functions do not fit models, access
# matrix storage, schedule workers, or mutate the S4 facade.
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

integrate_guides_current <- function(
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

assemble_pair_results <- function(raw_result, requested_pairs) {
    result <- data.table::rbindlist(
        list(raw_result, requested_pairs |> dplyr::filter(!pass_qc)), fill = TRUE
    )
    data.table::setorderv(result, c("p_value", "response_id"), na.last = TRUE)
    result
}

finalize_association_result <- function(result, integration_strategy, target_map,
                                        correction_method, alpha,
                                        is_calibration_check = FALSE) {
    result <- integrate_guides_current(result, integration_strategy, target_map,
                                       is_calibration_check)
    # Output amount 2/3 already has a skew-normal shape column named alpha.
    # Resolve correction settings outside the table's data mask.
    significant <- stats::p.adjust(result$p_value, method = correction_method) < alpha
    dplyr::mutate(result, significant = !!significant)
}

# Pinned pipeline scripts call these private names. They extract facade options
# and delegate; standard finalization follows gathering the complete family.
# Nuclear raw shards retain their separate existing script-level contract.
process_calibration_result <- function(result, sceptre_object) {
    finalize_association_result(result, sceptre_object@grna_integration_strategy,
        sceptre_object@grna_target_data_frame, sceptre_object@multiple_testing_method,
        sceptre_object@multiple_testing_alpha, TRUE)
}

process_discovery_result <- function(result, sceptre_object) {
    finalize_association_result(result, sceptre_object@grna_integration_strategy,
        sceptre_object@grna_target_data_frame, sceptre_object@multiple_testing_method,
        sceptre_object@multiple_testing_alpha)
}

apply_grouping_to_result <- function(result, sceptre_object,
                                     is_calibration_check = FALSE) {
    integrate_guides_current(result, sceptre_object@grna_integration_strategy,
        sceptre_object@grna_target_data_frame, is_calibration_check)
}

construct_data_frame_v2 <- function(curr_df, curr_response_result, output_amount) {
    pair_results_from_native(curr_df, curr_response_result, output_amount)
}
