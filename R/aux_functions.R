construct_data_frame_v2 <- function(curr_df, curr_response_result, return_debugging_metrics, return_resampling_dist) {
    curr_df$p_value <- sapply(X = curr_response_result, FUN = function(l) l$p, simplify = TRUE)
    curr_df$log_2_fold_change <- sapply(curr_response_result, FUN = function(l) l$lfc)
    if (return_debugging_metrics) {
      curr_df$stage <- sapply(curr_response_result, FUN = function(l) l$stage)
      curr_df$z_orig <- sapply(curr_response_result, FUN = function(l) l$z_orig)
      curr_df$xi <- sapply(curr_response_result, FUN = function(l) l$sn_params[1L])
      curr_df$omega <- sapply(curr_response_result, FUN = function(l) l$sn_params[2L])
      curr_df$alpha <- sapply(curr_response_result, FUN = function(l) l$sn_params[3L])
    }
    if (return_resampling_dist) {
      to_append <- lapply(curr_response_result, FUN = function(l) {
        m <- data.table::as.data.table(matrix(l$resampling_dist, nrow = 1))
        colnames(m) <- paste0("z_null_", seq(1L, ncol(m)))
        m
      }) |> data.table::rbindlist(fill = TRUE)
      curr_df <- cbind(curr_df, to_append)
    }
    return(curr_df)
}


#' Obtain discovery set
#'
#' The \code{obtain_discovery_set()} function returns the discovery set from an input discovery result.
#'
#' @inheritParams compare_calibration_and_discovery_results
#'
#' @return a subset of \code{discovery_result} containing the discoveries.
#' @export
#' @examples
#' # See the example in the run_sceptre_lowmoi help file.
#' ?run_sceptre_lowmoi
obtain_discovery_set <- function(discovery_result, alpha = 0.1, multiple_testing_correction = "BH") {
  discovery_result |>
    stats::na.omit() |>
    dplyr::mutate(p_adj = stats::p.adjust(p_value, method = multiple_testing_correction),
                  reject = p_adj < alpha) |>
    dplyr::filter(reject) |>
    dplyr::select(-p_adj, -reject) |>
    dplyr::arrange(p_value)
}


#' Generate all pairs
#'
#' The \code{generate_all_pairs()} function generates all possible response-gRNA group pairs that can be constructed from a given \code{response_matrix} and \code{grna_group_data_frame}.
#'
#' @inheritParams run_sceptre_lowmoi
#'
#' @return a data frame containing the columns \code{response_id} and \code{grna_group} in which each response is mapped to the entire set of targeting gRNA groups.
#' @export
#' @examples
#' data(response_matrix_lowmoi)
#' data(grna_group_data_frame_lowmoi)
#' response_grna_group_pairs <- generate_all_pairs(response_matrix_lowmoi,
#' grna_group_data_frame_lowmoi)
generate_all_pairs <- function(response_matrix, grna_group_data_frame) {
  response_ids <- rownames(response_matrix) |> factor()
  grna_groups <- grna_group_data_frame |>
    dplyr::filter(grna_group != "non-targeting") |>
    dplyr::pull(grna_group) |> unique() |> factor()
  expand.grid(response_id = response_ids, grna_group = grna_groups)
}
