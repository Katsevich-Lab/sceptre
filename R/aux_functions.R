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


auto_construct_formula_object <- function(cell_covariates, include_grna_covariates) {
  cell_covariate_names <- colnames(cell_covariates)
  cell_covariate_names <- cell_covariate_names[cell_covariate_names != "response_p_mito"]
  if (!include_grna_covariates) { # by default, do not use grna count-based covariates in low moi
    cell_covariate_names <- cell_covariate_names[!(cell_covariate_names %in% c("grna_n_umis", "grna_n_nonzero"))]
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
  form <- paste0("~ ", form_str) |> stats::as.formula()
  return(form)
}


auto_compute_cell_covariates <- function(response_matrix, grna_matrix, extra_covariates, response_names) {
  # compute the response covariates
  covariate_df <- compute_cell_covariates(response_matrix, response_names, TRUE)
  colnames(covariate_df) <- paste0("response_", colnames(covariate_df))

  # compute the grna covariates
  grna_covariate_df <- compute_cell_covariates(grna_matrix, response_names, FALSE)
  colnames(grna_covariate_df) <- paste0("grna_", colnames(grna_covariate_df))
  covariate_df <- cbind(covariate_df, grna_covariate_df)

  # if extra covariates have been provided, add those as well
  if (nrow(extra_covariates) >= 1L) {
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


get_log_dir <- function() {
  log_dir <- paste0(tempdir(), "/sceptre_logs/")
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
  return(log_dir)
}
