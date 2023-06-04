# helper function 0: compute D matrix
compute_D_matrix <- function(Zt_wZ, wZ) {
  P_decomp <- eigen(Zt_wZ, symmetric = TRUE)
  U <- P_decomp$vectors
  Lambda_minus_half <- 1/sqrt(P_decomp$values)
  D <- (Lambda_minus_half * t(U)) %*% t(wZ)
  return(D)
}

# helper function 1: idx vector for discovery analysis complement set
get_idx_vector_discovery_complement_set <- function(curr_grna_group, grna_group_idxs, n_cells) {
  trt_idxs <- grna_group_idxs[[curr_grna_group]]
  cntrl_idxs <- seq(1L, n_cells)[-trt_idxs]
  idxs <- c(cntrl_idxs, trt_idxs)
  return(list(idxs = idxs, n_trt = length(trt_idxs), n_cntrl = length(cntrl_idxs)))
}

# helper function 3: idx vector for calibration check complement set
get_idx_vector_calibration_check_complement_set <- function(curr_grna_group, indiv_nt_grna_idxs, n_cells, low_moi) {
  undercover_nts <- strsplit(x = curr_grna_group, split = "&", fixed = TRUE)[[1]]
  trt_idxs <- indiv_nt_grna_idxs[undercover_nts] |> unlist() |> stats::setNames(NULL)
  if (!low_moi) trt_idxs <- unique(trt_idxs)
  cntrl_idxs <- seq(1L, n_cells)[-trt_idxs]
  idxs <- c(cntrl_idxs, trt_idxs)
  return(list(idxs = idxs, n_trt = length(trt_idxs), n_cntrl = length(cntrl_idxs)))
}

# helper function 4: get undercover idx vector
get_idx_vector_calibration_check_nt_cells <- function(undercover_group, indiv_nt_grna_idxs) {
  undercover_nts <- strsplit(x = undercover_group, split = "&", fixed = TRUE)[[1]]
  control_cells <- indiv_nt_grna_idxs[setdiff(names(indiv_nt_grna_idxs), undercover_nts)] |>
    unlist() |> stats::setNames(NULL)
  undercover_cells <- indiv_nt_grna_idxs[undercover_nts] |> unlist() |> stats::setNames(NULL)
  idxs <- c(control_cells, undercover_cells)
  return(list(idxs = idxs, n_trt = length(undercover_cells), n_cntrl = length(control_cells)))
}

# workhorse function 2: discovery analysis nt cells perm test
discovery_ntcells_perm_test <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, covariate_matrix, all_nt_idxs, grna_group_idxs, indiv_nt_grna_idxs, grna_groups, pieces_precomp, expression_vector, response_precomp) {
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  n_cntrl <- length(all_nt_idxs)
  for (i in seq(grna_groups)) {
    grna_group <- grna_groups[i]
    trt_idxs <- grna_group_idxs[[grna_group]]
    n_trt <- length(trt_idxs)
    combined_idxs <- c(all_nt_idxs, trt_idxs)

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
    # 3. compute D matrix
    D <- compute_D_matrix(precomp_pieces$Zt_wZ, precomp_pieces$wZ)

    # 4. run the association test
    result <- run_low_level_test_full_v2(y = curr_expression_vector,
                                         mu = precomp_pieces$mu,
                                         a = precomp_pieces$a,
                                         w = precomp_pieces$w,
                                         D = D,
                                         n_cntrl = n_cntrl,
                                         n_trt = n_trt,
                                         synthetic_idxs = synthetic_idxs,
                                         B1 = B1, B2 = B2, B3 = B3,
                                         fit_skew_normal = fit_skew_normal,
                                         return_resampling_dist = return_resampling_dist)
    result_list_inner[[i]] <- result
  }
  return(result_list_inner)
}

# workhorse function 3: calibrabtion check complement set perm test
calibration_complement_perm_test <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, covariate_matrix, all_nt_idxs, grna_group_idxs, indiv_nt_grna_idxs, grna_groups, pieces_precomp, expression_vector, response_precomp) {
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  n_cells <- length(expression_vector)
  D <- compute_D_matrix(Zt_wZ = pieces_precomp$Zt_wZ, wZ = pieces_precomp$wZ)
  for (i in seq_along(grna_groups)) {
    curr_grna_group <- grna_groups[i]
    idxs <- get_idx_vector_calibration_check_complement_set(curr_grna_group = curr_grna_group,
                                                            indiv_nt_grna_idxs = indiv_nt_grna_idxs,
                                                            n_cells = n_cells, low_moi = TRUE)
    result <- run_low_level_test_full_v2(y = expression_vector[idxs$idxs],
                                         mu = pieces_precomp$mu[idxs$idxs],
                                         a = pieces_precomp$a[idxs$idxs],
                                         w = pieces_precomp$w[idxs$idxs],
                                         D = D[,idxs$idxs],
                                         n_cntrl = idxs$n_cntrl,
                                         n_trt = idxs$n_trt,
                                         synthetic_idxs = synthetic_idxs,
                                         B1 = B1, B2 = B2, B3 = B3,
                                         fit_skew_normal = fit_skew_normal,
                                         return_resampling_dist = return_resampling_dist)
    result_list_inner[[i]] <- result
  }
  return(result_list_inner)
}

# workhorse function 4: discovery analysis complement set perm test (different than the others)
discovery_complement_perm_test <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, covariate_matrix, all_nt_idxs, grna_group_idxs, indiv_nt_grna_idxs, grna_groups, pieces_precomp, expression_vector, response_precomp) {
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  n_cells <- length(expression_vector)
  D <- compute_D_matrix(Zt_wZ = pieces_precomp$Zt_wZ, wZ = pieces_precomp$wZ)
  for (i in seq_along(grna_groups)) {
    curr_grna_group <- grna_groups[i]
    idxs <- get_idx_vector_discovery_complement_set(curr_grna_group, grna_group_idxs, n_cells)
    result <- run_low_level_test_full_v2(y = expression_vector[idxs$idxs],
                                         mu = pieces_precomp$mu[idxs$idxs],
                                         a = pieces_precomp$a[idxs$idxs],
                                         w = pieces_precomp$w[idxs$idxs],
                                         D = D[,idxs$idxs],
                                         n_cntrl = idxs$n_cntrl,
                                         n_trt = idxs$n_trt,
                                         synthetic_idxs = synthetic_idxs,
                                         B1 = B1, B2 = B2, B3 = B3,
                                         fit_skew_normal = fit_skew_normal,
                                         return_resampling_dist = return_resampling_dist)
    result_list_inner[[i]] <- result
  }
  return(result_list_inner)
}
