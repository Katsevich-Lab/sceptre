highmoi_discovery_stat <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, grna_group_idxs, covariate_matrix, grna_groups, pieces_precomp, expression_vector) {
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  n_cells <- length(expression_vector)
  # compute the D matrix
  D <- compute_D_matrix(Zt_wZ = pieces_precomp$Zt_wZ, wZ = pieces_precomp$wZ)
  for (i in seq_along(grna_groups)) {
    curr_grna_group <- grna_groups[i]
    idxs <- get_discovery_idx_vector(curr_grna_group, grna_group_idxs, n_cells)
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


lowmoi_approximate_stat_discovery <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, grna_group_idxs, covariate_matrix, all_nt_idxs, regression_method, indiv_nt_grna_idxs, grna_groups, pieces_precomp, expression_vector_nt, expression_vector, response_precomp) {
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  for (i in seq_along(grna_groups)) {
    grna_group <- grna_groups[i]
    expression_vector_trt <- expression_vector[grna_group_idxs[[grna_group]]]
    covariate_matrix_trt <- covariate_matrix[grna_group_idxs[[grna_group]],]

    # 1. compute the precomputation pieces for the treatment cells
    pieces_trt <- compute_precomputation_pieces(expression_vector_trt, covariate_matrix_trt, response_precomp$fitted_coefs, response_precomp$theta, TRUE)

    # 2. compute the shared weighted covariate matrix
    Zt_wZ <- pieces_trt$Zt_wZ + pieces_precomp$Zt_wZ

    # 3. compute the D matrix
    D <- compute_D_matrix(Zt_wZ = Zt_wZ, wZ = rbind(pieces_precomp$wZ, pieces_trt$wZ))

    # 4. create the combined mu, y, and a vectors
    y <- c(expression_vector_nt, expression_vector_trt)
    mu <- c(pieces_precomp$mu, pieces_trt$mu)
    a <- c(pieces_precomp$a, pieces_trt$a)
    w <- c(pieces_precomp$w, pieces_trt$w)

    # 5. call low-level estimation and testing function
    n_cntrl <- length(expression_vector_nt)
    n_trt <- length(expression_vector_trt)
    result <- run_low_level_test_full_v2(y = y,
                                         mu = mu,
                                         a = a,
                                         w = w,
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


lowmoi_undercover_stat <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, grna_group_idxs, covariate_matrix, all_nt_idxs, regression_method, indiv_nt_grna_idxs, grna_groups, pieces_precomp, expression_vector_nt, expression_vector, response_precomp) {
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  # compute the D matrix
  D <- compute_D_matrix(Zt_wZ = pieces_precomp$Zt_wZ, wZ = pieces_precomp$wZ)
  for (i in seq_along(grna_groups)) {
    undercover_group <- grna_groups[i]
    idxs <- get_undercover_idx_vector(undercover_group, indiv_nt_grna_idxs)
    result <- run_low_level_test_full_v2(y = expression_vector_nt[idxs$idxs],
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


lowmoi_exact_stat_discovery <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, grna_group_idxs, covariate_matrix, all_nt_idxs, regression_method, indiv_nt_grna_idxs, grna_groups, pieces_precomp, expression_vector_nt, expression_vector, response_precomp) {
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  for (i in seq_along(grna_groups)) {
    grna_group <- grna_groups[i]
    trt_idxs <- grna_group_idxs[[grna_group]]
    n_trt <- length(trt_idxs)
    n_cntrl <- length(all_nt_idxs)
    combined_idxs <- c(all_nt_idxs, trt_idxs)

    # 1. get the expression vector and covariate matrix
    curr_expression_vector <- expression_vector[combined_idxs]
    curr_covariate_matrix <- covariate_matrix[combined_idxs,]

    # 2. perform the response precomputation
    response_precomp <- perform_response_precomputation(expressions = curr_expression_vector,
                                                        covariate_matrix = curr_covariate_matrix,
                                                        regression_method = regression_method)
    result_list_inner[[i]] <- tryCatch({
      # 2. get the precomp pieces
      precomp_pieces <- compute_precomputation_pieces(expression_vector = curr_expression_vector,
                                                      covariate_matrix = curr_covariate_matrix,
                                                      fitted_coefs = response_precomp$fitted_coefs,
                                                      theta = response_precomp$theta,
                                                      full_test_stat = TRUE)
      # 3. compute D matrix
      D <- compute_D_matrix(precomp_pieces$Zt_wZ, precomp_pieces$wZ)

      # 4. run the association test
      run_low_level_test_full_v2(y = curr_expression_vector,
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
    }, error = function(e) backup_distilled(curr_expression_vector, curr_covariate_matrix,
                                            response_precomp, n_cntrl, n_trt, synthetic_idxs,
                                            B1, B2, B3, fit_skew_normal, return_resampling_dist),
    warning = function(w) backup_distilled(curr_expression_vector, curr_covariate_matrix,
                                           response_precomp, n_cntrl, n_trt, synthetic_idxs,
                                           B1, B2, B3, fit_skew_normal, return_resampling_dist))
  }
  return(result_list_inner)
}


backup_distilled <- function(curr_expression_vector, curr_covariate_matrix,
                             response_precomp, n_cntrl, n_trt, synthetic_idxs,
                             B1, B2, B3, fit_skew_normal, return_resampling_dist) {
  # 2. get the precomp pieces
  precomp_pieces <- compute_precomputation_pieces(expression_vector = curr_expression_vector,
                                                  covariate_matrix = curr_covariate_matrix,
                                                  fitted_coefs = response_precomp$fitted_coefs,
                                                  theta = response_precomp$theta,
                                                  full_test_stat = FALSE)
  # 4. run the association test
  run_low_level_test_distilled(y = curr_expression_vector,
                               mu = precomp_pieces$mu,
                               a = precomp_pieces$a,
                               b = precomp_pieces$b,
                               n_cntrl = n_cntrl,
                               n_trt = n_trt,
                               synthetic_idxs = synthetic_idxs,
                               B1 = B1, B2 = B2, B3 = B3,
                               fit_skew_normal = fit_skew_normal,
                               return_resampling_dist = return_resampling_dist)
}


compute_D_matrix <- function(Zt_wZ, wZ) {
  P_decomp <- eigen(Zt_wZ, symmetric = TRUE)
  U <- P_decomp$vectors
  Lambda_minus_half <- 1/sqrt(P_decomp$values)
  D <- (Lambda_minus_half * t(U)) %*% t(wZ)
  return(D)
}


get_undercover_idx_vector <- function(undercover_group, indiv_nt_grna_idxs) {
  undercover_nts <- strsplit(x = undercover_group, split = "&", fixed = TRUE)[[1]]
  control_cells <- indiv_nt_grna_idxs[setdiff(names(indiv_nt_grna_idxs), undercover_nts)] |>
    unlist() |> stats::setNames(NULL)
  undercover_cells <- indiv_nt_grna_idxs[undercover_nts] |> unlist() |> stats::setNames(NULL)
  idxs <- c(control_cells, undercover_cells)
  return(list(idxs = idxs, n_trt = length(undercover_cells), n_cntrl = length(control_cells)))
}


get_discovery_idx_vector <- function(curr_grna_group, grna_group_idxs, n_cells) {
  trt_idxs <- grna_group_idxs[[curr_grna_group]]
  cntrl_idxs <- seq(1L, n_cells)[-trt_idxs]
  idxs <- c(cntrl_idxs, trt_idxs)
  return(list(idxs = idxs, n_trt = length(trt_idxs), n_cntrl = length(cntrl_idxs)))
}
