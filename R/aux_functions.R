construct_data_frame_v2 <- function(curr_df, curr_response_result, output_amount) {
    curr_df$p_value <- sapply(X = curr_response_result, FUN = function(l) l$p, simplify = TRUE)
    curr_df$log_2_fold_change <- sapply(curr_response_result, FUN = function(l) l$lfc)
    if (output_amount >= 2L) {
      curr_df$stage <- sapply(curr_response_result, FUN = function(l) l$stage)
      curr_df$z_orig <- sapply(curr_response_result, FUN = function(l) l$z_orig)
      curr_df$xi <- sapply(curr_response_result, FUN = function(l) l$sn_params[1L])
      curr_df$omega <- sapply(curr_response_result, FUN = function(l) l$sn_params[2L])
      curr_df$alpha <- sapply(curr_response_result, FUN = function(l) l$sn_params[3L])
    }
    if (output_amount >= 3L) {
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
generate_all_pairs <- function(response_matrix, grna_group_data_frame) {
  response_ids <- rownames(response_matrix) |> factor()
  grna_groups <- grna_group_data_frame |>
    dplyr::filter(grna_group != "non-targeting") |>
    dplyr::pull(grna_group) |> unique() |> factor()
  expand.grid(response_id = response_ids, grna_group = grna_groups)
}


auto_construct_formula_object <- function(cell_covariates, low_moi) {
  cell_covariate_names <- colnames(cell_covariates)
  if (low_moi) { # by default, do not use grna count-based covariates in low moi
    cell_covariate_names <- cell_covariate_names[cell_covariate_names != c("grna_n_umis", "grna_n_nonzero")]
  }
  form_str <- sapply(cell_covariate_names, function(curr_name) {
    count_based_covariate <- grepl(pattern = "n_umis|n_nonzero", x = curr_name)
    if (count_based_covariate) {
      if (any(cell_covariates[[curr_name]] == 0)) {
        out <- paste0("log(", curr_name, "+1)")
      } else {
        out <- paste0("log(", curr_name, ")")
      }
    } else {
      out <- curr_name
    }
    return(out)
  }) |> paste0(collapse = " + ")
  form <- paste0("~ ", form_str) |> as.formula()
  return(form)
}


auto_compute_cell_covariates <- function(response_matrix, grna_matrix, extra_covariates) {
  # compute the response covariates
  covariate_df <- compute_cell_covariates(response_matrix)
  colnames(covariate_df) <- paste0("response_", colnames(covariate_df))

  # compute the grna covariates
  grna_covariate_df <- compute_cell_covariates(grna_matrix)
  colnames(grna_covariate_df) <- paste0("grna_", colnames(grna_covariate_df))
  covariate_df <- cbind(covariate_df, grna_covariate_df)

  # if extra covariates have been provided, add those as well
  if (!is.null(extra_covariates)) {
    covariate_df <- cbind(covariate_df, extra_covariates)
  }
  return(covariate_df)
}


partition_response_ids <- function(response_ids, parallel) {
  groups_set <- FALSE
  if (parallel) {
    n_cores <- parallel::detectCores(logical = FALSE)
    if (length(response_ids) >= n_cores) {
      set.seed(4)
      s <- sample(response_ids)
      out <- split(s, cut(seq_along(s), n_cores, labels = paste0("group_", seq(1, n_cores))))
      groups_set <- TRUE
    }
  }
  if (!groups_set) out <- list(group_1 = response_ids)
  return(out)
}
