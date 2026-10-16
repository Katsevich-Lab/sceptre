# Ordinary raw-table operations; these functions do not fit models, access
# matrix storage, schedule workers, or mutate the S4 facade.
#' Attach aligned native records to requested pair rows
#'
#' One native record per row supplies p/fold-change/SE. Output amount >=2 adds
#' stage/score/tail parameters; >=3 binds full null draws with fill for unequal
#' lengths. Retains pair rows and the existing data-frame/table binding classes.
#' No guide integration or multiple-testing correction occurs at this boundary.
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

#' Integrate guide rows using the configured existing result rule
#'
#' Union relabels groups as targets. Singleton/Bonferroni join the guide-target
#' map, retaining its multiplicity, then sort by p-value/response. Non-calibration
#' Bonferroni retains one minimum-p record per response/target and multiplies by
#' passing-guide count, capped at one; all-failed groups keep maximum counts and
#' missing p/effect. Calibration skips that collapse. No new correction occurs.
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

#' Restore failed-QC requests before final family processing
#'
#' Binds tested rows followed by failed requests, filling missing native fields,
#' then orders by p-value/response with NA last. Returns a data.table. Requested
#' identities and repeated rows remain; no p-values are computed or adjusted.
#' @noRd
assemble_pair_results <- function(raw_result, requested_pairs) {
    result <- data.table::rbindlist(
        list(raw_result, requested_pairs |> dplyr::filter(!pass_qc)), fill = TRUE
    )
    data.table::setorderv(result, c("p_value", "response_id"), na.last = TRUE)
    result
}

#' Integrate the complete result family and assign its existing significance flag
#'
#' Uses p.adjust with the facade's method and strict adjusted-p < alpha after
#' guide integration. Resolve alpha outside the table mask: amount 2/3 carries
#' a separate skew-normal shape column with the same name. Real pipeline callers
#' invoke this after gathering the full family, rather than per raw shard.
#' @noRd
finalize_association_result <- function(result, integration_strategy, target_map,
                                        correction_method, alpha,
                                        is_calibration_check = FALSE) {
    result <- integrate_guides(result, integration_strategy, target_map,
                                       is_calibration_check)
    # Output amount 2/3 already has a skew-normal shape column named alpha.
    # Resolve correction settings outside the table's data mask.
    significant <- stats::p.adjust(result$p_value, method = correction_method) < alpha
    dplyr::mutate(result, significant = !!significant)
}

# Pinned pipeline scripts call these private names. They extract facade options
# and delegate; standard finalization follows gathering the complete family.
# Nuclear raw shards retain their separate existing script-level contract.
#' Finalize a gathered calibration family using facade settings
#'
#' Retains the direct Nextflow private-name contract and skips the non-calibration
#' Bonferroni collapse used by power/discovery. Returns the finalized table;
#' does not update facade slots.
#' @noRd
process_calibration_result <- function(result, sceptre_object) {
    finalize_association_result(result, sceptre_object@grna_integration_strategy,
        sceptre_object@grna_target_data_frame, sceptre_object@multiple_testing_method,
        sceptre_object@multiple_testing_alpha, TRUE)
}

#' Finalize a gathered discovery family using facade settings
#'
#' Direct Nextflow adapter for the complete-family correction above; returns a
#' table without changing facade slots or recomputing tests.
#' @noRd
process_discovery_result <- function(result, sceptre_object) {
    finalize_association_result(result, sceptre_object@grna_integration_strategy,
        sceptre_object@grna_target_data_frame, sceptre_object@multiple_testing_method,
        sceptre_object@multiple_testing_alpha)
}

#' Apply guide grouping alone at the actual pipeline boundary
#'
#' Returns integrated rows using facade metadata and the explicit calibration
#' flag. This adapter deliberately does not assign complete-family significance;
#' existing nuclear scripts retain their separate raw-shard processing contract.
#' @noRd
apply_grouping_to_result <- function(result, sceptre_object,
                                     is_calibration_check = FALSE) {
    integrate_guides(result, sceptre_object@grna_integration_strategy,
        sceptre_object@grna_target_data_frame, is_calibration_check)
}
