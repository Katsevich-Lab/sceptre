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


lowmoi_distilled_stat_undercover <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, indiv_nt_grna_idxs, grna_groups, pieces_precomp, expression_vector_nt) {
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  for (i in seq_along(grna_groups)) {
    undercover_group <- grna_groups[i]
    idxs <- get_undercover_idx_vector(undercover_group, indiv_nt_grna_idxs)
    result <- run_low_level_test_distilled(y = expression_vector_nt[idxs$idxs],
                                           mu = pieces_precomp$mu[idxs$idxs],
                                           a = pieces_precomp$a[idxs$idxs],
                                           b = pieces_precomp$b[idxs$idxs],
                                           n_cntrl = idxs$n_cntrl, n_trt = idxs$n_trt,
                                           synthetic_idxs = synthetic_idxs,
                                           B1 = B1, B2 = B2, B3 = B3,
                                           fit_skew_normal = fit_skew_normal,
                                           return_resampling_dist = return_resampling_dist)
    result_list_inner[[i]] <- result
  }
  return(result_list_inner)
}


lowmoi_full_stat_undercover <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, indiv_nt_grna_idxs, grna_groups, pieces_precomp, expression_vector_nt) {
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


lowmoi_full_stat_discovery <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, grna_group_idxs, covariate_matrix, grna_groups, pieces_precomp, expression_vector_nt, expression_vector, response_precomp) {
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


lowmoi_distilled_stat_discovery <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, grna_group_idxs, covariate_matrix, grna_groups, pieces_precomp, expression_vector_nt, expression_vector, response_precomp) {
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  for (i in seq_along(grna_groups)) {
    grna_group <- grna_groups[i]
    expression_vector_trt <- expression_vector[grna_group_idxs[[grna_group]]]
    covariate_matrix_trt <- covariate_matrix[grna_group_idxs[[grna_group]],]
    
    # 1. compute the pieces
    pieces_trt <- compute_precomputation_pieces(expression_vector_trt, covariate_matrix_trt, response_precomp$fitted_coefs, response_precomp$theta, FALSE)
    
    # 2. set the arguments
    y <- c(expression_vector_nt, expression_vector_trt)
    mu <- c(pieces_precomp$mu, pieces_trt$mu)
    a <- c(pieces_precomp$a, pieces_trt$a)
    b <- c(pieces_precomp$b, pieces_trt$b)
    
    # 3. call function
    n_cntrl <- length(expression_vector_nt)
    n_trt <- length(expression_vector_trt)
    result <- run_low_level_test_distilled(y = y, mu = mu, a = a, b = b,
                                           n_cntrl = n_cntrl, n_trt = n_trt,
                                           synthetic_idxs = synthetic_idxs,
                                           B1 = B1, B2 = B2, B3 = B3,
                                           fit_skew_normal = fit_skew_normal,
                                           return_resampling_dist = return_resampling_dist)
    result_list_inner[[i]] <- result
  }
  return(result_list_inner)
}
