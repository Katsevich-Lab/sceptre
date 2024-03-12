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

# helper function 1. calibration check complement set
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


# workhorse function 1. permutations, glm factored out
perm_test_glm_factored_out <- function(synthetic_idxs, B1, B2, B3, fit_parametric_curve, output_amount, grna_groups,
                                       expression_vector, pieces_precomp, get_idx_f, side_code) {
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
                                         fit_parametric_curve = fit_parametric_curve,
                                         return_resampling_dist = (output_amount == 3L),
                                         side_code = side_code)
    result_list_inner[[i]] <- result
  }
  return(result_list_inner)
}


# workhorse function 2: permutations, glm run inside
discovery_ntcells_perm_test <- function(synthetic_idxs, B1, B2, B3, fit_parametric_curve, output_amount, covariate_matrix,
                                        all_nt_idxs, grna_group_idxs, grna_groups, expression_vector, side_code) {
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  for (i in seq_along(grna_groups)) {
    curr_grna_group <- grna_groups[i]
    idxs <- get_idx_vector_discovery_analysis(curr_grna_group = curr_grna_group, grna_group_idxs = grna_group_idxs)
    n_trt <- idxs$n_trt
    combined_idxs <- c(idxs$trt_idxs, all_nt_idxs)
    trt_idxs <- seq(1L, n_trt)

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
                                         trt_idxs = trt_idxs,
                                         synthetic_idxs = synthetic_idxs,
                                         B1 = B1, B2 = B2, B3 = B3,
                                         fit_parametric_curve = fit_parametric_curve,
                                         return_resampling_dist = (output_amount == 3L),
                                         side_code = side_code)
    result_list_inner[[i]] <- result
  }
  return(result_list_inner)
}


# workhorse function 3: crt, glm factored out
crt_glm_factored_out <- function(B1, B2, B3, fit_parametric_curve, output_amount,
                                 response_ids, response_precomputations, covariate_matrix,
                                 get_idx_f, curr_grna_group, subset_to_nt_cells, all_nt_idxs,
                                 response_matrix, side_code, cells_in_use) {
  result_list_inner <- vector(mode = "list", length = length(response_ids))
  # precomputation on grna
  idxs <- get_idx_f(curr_grna_group)
  trt_idxs <- idxs$trt_idxs
  n_trt <- idxs$n_trt
  fitted_probabilities <- perform_grna_precomputation(trt_idxs = trt_idxs,
                                                      covariate_matrix = covariate_matrix,
                                                      return_fitted_values = TRUE)
  synthetic_idxs <- crt_index_sampler_fast(fitted_probabilities = fitted_probabilities, B = B1 + B2 + B3)
  # loop over genes
  for (i in seq_along(response_ids)) {
    curr_response_id <- response_ids[i]
    # load gene expressions
    expression_vector <- load_csr_row(j = response_matrix@j,
                                      p = response_matrix@p,
                                      x = response_matrix@x,
                                      row_idx = which(rownames(response_matrix) == curr_response_id),
                                      n_cells = ncol(response_matrix))[cells_in_use]
    if (subset_to_nt_cells) expression_vector <- expression_vector[all_nt_idxs]
    # compute the precomputation pieces
    pieces_precomp <- compute_precomputation_pieces(expression_vector = expression_vector,
                                                    covariate_matrix = covariate_matrix,
                                                    fitted_coefs = response_precomputations[[curr_response_id]]$fitted_coefs,
                                                    theta = response_precomputations[[curr_response_id]]$theta,
                                                    full_test_stat = TRUE)
    # run the association test
    result <- run_low_level_test_full_v4(y = expression_vector,
                                         mu = pieces_precomp$mu,
                                         a = pieces_precomp$a,
                                         w = pieces_precomp$w,
                                         D = pieces_precomp$D,
                                         trt_idxs = trt_idxs,
                                         n_trt = n_trt,
                                         use_all_cells = TRUE,
                                         synthetic_idxs = synthetic_idxs,
                                         B1 = B1, B2 = B2, B3 = B3,
                                         fit_parametric_curve = fit_parametric_curve,
                                         return_resampling_dist = (output_amount == 3L),
                                         side_code = side_code)
    result_list_inner[[i]] <- result
  }
  return(result_list_inner)
}

# workhorse function 4: crt, glm run inside
discovery_ntcells_crt <- function(B1, B2, B3, fit_parametric_curve, output_amount, get_idx_f, response_ids,
                                  covariate_matrix, curr_grna_group, all_nt_idxs, response_matrix,
                                  side_code, cells_in_use) {
  result_list_inner <- vector(mode = "list", length = length(response_ids))
  # initialize the idxs
  idxs <- get_idx_f(curr_grna_group)
  n_trt <- idxs$n_trt
  combined_idxs <- c(idxs$trt_idxs, all_nt_idxs)
  trt_idxs <- seq(1L, n_trt)
  # subset the covariate matrix for this gRNA
  curr_covariate_matrix <- covariate_matrix[combined_idxs,]
  # perform the grna precomputation
  fitted_probabilities <- perform_grna_precomputation(trt_idxs = trt_idxs,
                                                      covariate_matrix = curr_covariate_matrix,
                                                      return_fitted_values = TRUE)
  synthetic_idxs <- crt_index_sampler_fast(fitted_probabilities = fitted_probabilities, B = B1 + B2 + B3)
  # loop over the response ids
  for (i in seq_along(response_ids)) {
    curr_response_id <- response_ids[i]
    # load gene expressions
    curr_expression_vector <- load_csr_row(j = response_matrix@j,
                                           p = response_matrix@p,
                                           x = response_matrix@x,
                                           row_idx = which(rownames(response_matrix) == curr_response_id),
                                           n_cells = ncol(response_matrix))[cells_in_use]
    curr_expression_vector <- curr_expression_vector[combined_idxs]
    # perform the response precomputation
    response_precomp <- perform_response_precomputation(expressions = curr_expression_vector,
                                                        covariate_matrix = curr_covariate_matrix,
                                                        regression_method = "pois_glm")
    # get the precomp pieces
    pieces_precomp <- compute_precomputation_pieces(expression_vector = curr_expression_vector,
                                                    covariate_matrix = curr_covariate_matrix,
                                                    fitted_coefs = response_precomp$fitted_coefs,
                                                    theta = response_precomp$theta,
                                                    full_test_stat = TRUE)
    # run the association test
    result <- run_low_level_test_full_v4(y = curr_expression_vector,
                                         mu = pieces_precomp$mu,
                                         a = pieces_precomp$a,
                                         w = pieces_precomp$w,
                                         D = pieces_precomp$D,
                                         trt_idxs = trt_idxs,
                                         n_trt = n_trt,
                                         use_all_cells = TRUE,
                                         synthetic_idxs = synthetic_idxs,
                                         B1 = B1, B2 = B2, B3 = B3,
                                         fit_parametric_curve = fit_parametric_curve,
                                         return_resampling_dist = (output_amount == 3L),
                                         side_code = side_code)
    result_list_inner[[i]] <- result
  }
  return(result_list_inner)
}

# workhorse function 5: score test, glm factored out
score_test_glm_factored_out <- function(grna_groups, expression_vector,
                                        pieces_precomp, get_idx_f, side_code){
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  for (i in seq_along(grna_groups)) {
    curr_grna_group <- grna_groups[i]
    idxs <- get_idx_f(curr_grna_group)

    z_orig = compute_observed_full_statistic_v2(a = pieces_precomp$a,
                                                w = pieces_precomp$w,
                                                D = pieces_precomp$D,
                                                trt_idxs = idxs$trt_idxs)

    result_list_inner[[i]] <- result
  }
  return(result_list_inner)
}
