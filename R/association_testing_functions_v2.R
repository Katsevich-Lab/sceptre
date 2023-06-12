# helper function 0: get_idx_vector function factory
get_idx_vector_factory <- function(calibration_check, indiv_nt_grna_idxs, grna_group_idxs, low_moi) {
    f <- function(curr_grna_group) {
      if (calibration_check) {
        get_idx_vector_calibration_check(curr_grna_group, indiv_nt_grna_idxs, low_moi)
      } else {
        get_idx_vector_discovery_analysis(curr_grna_group, grna_group_idxs)
      }
    }
    return(f)
}

# helper function 1: calibration check complement set
get_idx_vector_calibration_check <- function(curr_grna_group, indiv_nt_grna_idxs, low_moi) {
  undercover_nts <- strsplit(x = curr_grna_group, split = "&", fixed = TRUE)[[1]]
  trt_idxs <- indiv_nt_grna_idxs[undercover_nts] |> unlist() |> stats::setNames(NULL)
  if (!low_moi) trt_idxs <- unique(trt_idxs)
  return(list(trt_idxs = trt_idxs, n_trt = length(trt_idxs)))
}

# helper function 2: discovery analysis
get_idx_vector_discovery_analysis <- function(curr_grna_group, grna_group_idxs) {
  trt_idxs <- grna_group_idxs[[curr_grna_group]]
  return(list(trt_idxs = trt_idxs, n_trt = length(trt_idxs)))
}

# permutation workhorse function 1 (glm factored out)
perm_test_glm_factored_out <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, grna_groups, expression_vector, pieces_precomp, get_idx_f) {
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  for (i in seq_along(grna_groups)) {
    curr_grna_group <- grna_groups[i]
    idxs <- get_idx_f(curr_grna_group)
    result <- run_low_level_test_full_v4(y = expression_vector,
                                         mu = pieces_precomp$mu,
                                         a = pieces_precomp$a,
                                         w = pieces_precomp$w,
                                         D = pieces_precomp$D,
                                         trt_idxs = idxs$trt_idxs,
                                         use_all_cells = FALSE,
                                         n_trt = idxs$n_trt,
                                         synthetic_idxs = synthetic_idxs,
                                         B1 = B1, B2 = B2, B3 = B3,
                                         fit_skew_normal = fit_skew_normal,
                                         return_resampling_dist = return_resampling_dist)
    result_list_inner[[i]] <- result
  }
  return(result_list_inner)
}


# workhorse function 2 (glm to be run)
discovery_ntcells_perm_test <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, covariate_matrix, all_nt_idxs, grna_group_idxs, grna_groups, expression_vector) {
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  for (i in seq_along(grna_groups)) {
    curr_grna_group <- grna_groups[i]
    idxs <- get_idx_vector_discovery_analysis(curr_grna_group = curr_grna_group, grna_group_idxs = grna_group_idxs)
    trt_idxs <- idxs$trt_idxs
    n_trt <- idxs$n_trt
    combined_idxs <- c(trt_idxs, all_nt_idxs)

    # 1. get the expression vector and covariate matrix
    curr_expression_vector <- expression_vector[combined_idxs]
    curr_covariate_matrix <- covariate_matrix[combined_idxs,]

    # 2. perform the response precomputation
    response_precomp <- perform_response_precomputation(expressions = curr_expression_vector,
                                                        covariate_matrix = curr_covariate_matrix,
                                                        regression_method = "pois_glm")
    # 3. get the precomp pieces
    precomp_pieces <- compute_precomputation_pieces(expression_vector = curr_expression_vector,
                                                    covariate_matrix = curr_covariate_matrix,
                                                    fitted_coefs = response_precomp$fitted_coefs,
                                                    theta = response_precomp$theta,
                                                    full_test_stat = TRUE)
    # 4. run the association test
    result <- run_low_level_test_full_v4(y = curr_expression_vector,
                                         mu = precomp_pieces$mu,
                                         a = precomp_pieces$a,
                                         w = precomp_pieces$w,
                                         D = precomp_pieces$D,
                                         n_trt = n_trt,
                                         use_all_cells = FALSE,
                                         trt_idxs = seq(1L, n_trt),
                                         synthetic_idxs = synthetic_idxs,
                                         B1 = B1, B2 = B2, B3 = B3,
                                         fit_skew_normal = fit_skew_normal,
                                         return_resampling_dist = return_resampling_dist)
    result_list_inner[[i]] <- result
  }
  return(result_list_inner)
}
