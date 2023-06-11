# helper function 1: calibration check complement set
get_idx_vector_calibration_check <- function(curr_grna_group, indiv_nt_grna_idxs, take_unique) {
  undercover_nts <- strsplit(x = curr_grna_group, split = "&", fixed = TRUE)[[1]]
  trt_idxs <- indiv_nt_grna_idxs[undercover_nts] |> unlist() |> stats::setNames(NULL)
  if (take_unique) trt_idxs <- unique(trt_idxs)
  return(list(trt_idxs = trt_idxs, n_trt = length(trt_idxs)))
}

# helper function 2: discovery analysis
get_idx_vector_discovery_analysis <- function(curr_grna_group, grna_group_idxs) {
  trt_idxs <- grna_group_idxs[[curr_grna_group]]
  return(list(trt_idxs = trt_idxs, n_trt = length(trt_idxs)))
}

# workhorse function 1: calibrabtion check complement set perm test
calibration_complement_perm_test <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, covariate_matrix, all_nt_idxs, grna_group_idxs, indiv_nt_grna_idxs, grna_groups, expression_vector, pieces_precomp) {
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  for (i in seq_along(grna_groups)) {
    curr_grna_group <- grna_groups[i]
    idxs <- get_idx_vector_calibration_check(curr_grna_group = curr_grna_group,
                                             indiv_nt_grna_idxs = indiv_nt_grna_idxs,
                                             take_unique = FALSE)
    result <- run_low_level_test_full_v3(y = expression_vector,
                                         mu = pieces_precomp$mu,
                                         a = pieces_precomp$a,
                                         w = pieces_precomp$w,
                                         D = pieces_precomp$D,
                                         trt_idxs = idxs$trt_idxs,
                                         n_trt = idxs$n_trt,
                                         synthetic_idxs = synthetic_idxs,
                                         B1 = B1, B2 = B2, B3 = B3,
                                         fit_skew_normal = fit_skew_normal,
                                         return_resampling_dist = return_resampling_dist)
    result_list_inner[[i]] <- result
  }
  return(result_list_inner)
}

# workhorse function 2: calibration_ntcells_perm_test
calibration_ntcells_perm_test <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, covariate_matrix, all_nt_idxs, grna_group_idxs, indiv_nt_grna_idxs, grna_groups, expression_vector, pieces_precomp) {
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  for (i in seq_along(grna_groups)) {
    curr_grna_group <- grna_groups[i]
    idxs <- get_idx_vector_calibration_check(curr_grna_group, indiv_nt_grna_idxs, take_unique = FALSE)
    result <- run_low_level_test_full_v3(y = expression_vector,
                                         mu = pieces_precomp$mu,
                                         a = pieces_precomp$a,
                                         w = pieces_precomp$w,
                                         D = pieces_precomp$D,
                                         n_trt = idxs$n_trt,
                                         trt_idxs = idxs$trt_idxs,
                                         synthetic_idxs = synthetic_idxs,
                                         B1 = B1, B2 = B2, B3 = B3,
                                         fit_skew_normal = fit_skew_normal,
                                         return_resampling_dist = return_resampling_dist)
    result_list_inner[[i]] <- result
  }
  return(result_list_inner)
}

# workhorse function 3: discovery analysis complement set perm test
discovery_complement_perm_test <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, covariate_matrix, all_nt_idxs, grna_group_idxs, indiv_nt_grna_idxs, grna_groups, expression_vector, pieces_precomp) {
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  for (i in seq_along(grna_groups)) {
    curr_grna_group <- grna_groups[i]
    idxs <- get_idx_vector_discovery_analysis(curr_grna_group, grna_group_idxs)
    result <- run_low_level_test_full_v3(y = expression_vector,
                                         mu = pieces_precomp$mu,
                                         a = pieces_precomp$a,
                                         w = pieces_precomp$w,
                                         D = pieces_precomp$D,
                                         trt_idxs = idxs$trt_idxs,
                                         n_trt = idxs$n_trt,
                                         synthetic_idxs = synthetic_idxs,
                                         B1 = B1, B2 = B2, B3 = B3,
                                         fit_skew_normal = fit_skew_normal,
                                         return_resampling_dist = return_resampling_dist)
    result_list_inner[[i]] <- result
  }
  return(result_list_inner)
}

# workhorse function 4: discovery analysis nt cells perm test
discovery_ntcells_perm_test <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, covariate_matrix, all_nt_idxs, grna_group_idxs, indiv_nt_grna_idxs, grna_groups, expression_vector, pieces_precomp) {
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
    result <- run_low_level_test_full_v3(y = curr_expression_vector,
                                         mu = precomp_pieces$mu,
                                         a = precomp_pieces$a,
                                         w = precomp_pieces$w,
                                         D = precomp_pieces$D,
                                         n_trt = n_trt,
                                         trt_idxs = seq(1L, n_trt),
                                         synthetic_idxs = synthetic_idxs,
                                         B1 = B1, B2 = B2, B3 = B3,
                                         fit_skew_normal = fit_skew_normal,
                                         return_resampling_dist = return_resampling_dist)
    result_list_inner[[i]] <- result
  }
  return(result_list_inner)
}
