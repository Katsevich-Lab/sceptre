# helper function 1
get_id_from_idx <- function(response_idx, print_progress, response_ids, print_multiple = 5L, gc_multiple = 200L) {
  if ((response_idx == 1 || response_idx %% print_multiple == 0) && print_progress) {
    cat(paste0("Analyzing pairs containing response ", as.character(response_ids[response_idx]),
               " (", response_idx, " of ", length(response_ids), ")\n"))
  }
  if (response_idx %% gc_multiple == 0) gc() |> invisible()
  response_id <- as.character(response_ids[response_idx])
  return(response_id)
}

# core function 1: run permutation test in memory
run_perm_test_in_memory <- function(response_matrix, grna_assignments,
                                    covariate_matrix, response_grna_group_pairs,
                                    synthetic_idxs, return_resampling_dist, fit_skew_normal,
                                    B1, B2, B3, calibration_check, control_group, n_nonzero_trt_thresh,
                                    n_nonzero_cntrl_thresh, return_debugging_metrics, print_progress) {
  # 0. define several variables
  result_list_outer <- vector(mode = "list", length = 2 * length(unique(response_grna_group_pairs$response_id)))
  out_counter <- 1L
  subset_to_nt_cells <- calibration_check && !control_group_complement
  run_outer_regression <- calibration_check || control_group_complement
  all_nt_idxs <- if (!is.null(grna_assignments$all_nt_idxs)) grna_assignments$all_nt_idxs else NA
  grna_group_idxs <- if (!is.null(grna_assignments$grna_group_idxs)) grna_assignments$grna_group_idxs else NA
  indiv_nt_grna_idxs <- if (!is.null(grna_assignments$indiv_nt_grna_idxs)) grna_assignments$indiv_nt_grna_idxs else NA
  n_cells <- nrow(covariate_matrix)
  get_idx_f <- get_idx_vector_factory(calibration_check, indiv_nt_grna_idxs, grna_group_idxs, low_moi)

  # 1. subset covariate matrix to nt cells (if applicable)
  if (subset_to_nt_cells) covariate_matrix <- covariate_matrix[all_nt_idxs,]

  # 2. loop over the response IDs
  response_ids <- unique(response_grna_group_pairs$response_id)
  for (response_idx in seq_along(response_ids)) {
    response_id <- get_id_from_idx(response_idx, print_progress, response_ids)

    # 3. load the expressions of the current response
    expression_vector <- load_csr_row(j = response_matrix@j,
                                      p = response_matrix@p,
                                      x = response_matrix@x,
                                      row_idx = which(rownames(response_matrix) == response_id),
                                      n_cells = n_cells)
    if (low_level_association_funct == "calibration_ntcells_perm_test") {
      expression_vector <- expression_vector[all_nt_idxs]
    }

    # 4. obtain the gRNA groups to analyze
    l <- response_grna_group_pairs$response_id == response_id
    curr_df <- response_grna_group_pairs[l,]

    # 5. run QC if conducting a discovery analysis
    if (!calibration_check) {
      gene_wise_qc_result <- do_genewise_qc(expression_vector, all_nt_idxs, control_group_complement, curr_df, grna_group_idxs, n_nonzero_trt_thresh, n_nonzero_cntrl_thresh)
      curr_df <- gene_wise_qc_result$curr_df
      pass_qc <- gene_wise_qc_result$pass_qc
      any_pass_qc <- gene_wise_qc_result$any_pass_qc
      # if no row passes qc, then jump to next iteration
      if (!any_pass_qc) {
        result_list_outer[[out_counter]] <- curr_df
        out_counter <- out_counter + 1L
        next
      } else {
        # if some rows pass qc, remove those that fail to pass qc, and keep those that do
        result_list_outer[[out_counter]] <- curr_df[!pass_qc,]
        out_counter <- out_counter + 1L
      }
      curr_df <- curr_df[pass_qc,]
    }
    grna_groups <- as.character(curr_df$grna_group)

    # 6. perform outer regression (if applicable)
    if (run_outer_regression) {
      response_precomp <- perform_response_precomputation(expressions = expression_vector,
                                                          covariate_matrix = covariate_matrix,
                                                          regression_method = "pois_glm")
      pieces_precomp <- compute_precomputation_pieces(expression_vector,
                                                      covariate_matrix,
                                                      response_precomp$fitted_coefs,
                                                      response_precomp$theta,
                                                      full_test_stat = TRUE)
      curr_response_result <- perm_test_glm_factored_out(synthetic_idxs, B1, B2, B3, fit_skew_normal,
                                                         return_resampling_dist, grna_groups,
                                                         expression_vector, pieces_precomp, get_idx_f)
    } else {
      curr_response_result <- discovery_ntcells_perm_test(synthetic_idxs, B1, B2, B3, fit_skew_normal,
                                                          return_resampling_dist, covariate_matrix, all_nt_idxs,
                                                          grna_group_idxs, grna_groups, expression_vector)
    }

    # 9. combine the response-wise results into a data table; insert into list
    result_list_outer[[out_counter]] <- construct_data_frame_v2(curr_df, curr_response_result,
                                                                return_debugging_metrics, return_resampling_dist)
    out_counter <- out_counter + 1L
  }

  # combine and sort result
  ret <- data.table::rbindlist(result_list_outer, fill = TRUE)
  data.table::setorderv(ret, cols = c("p_value", "response_id"), na.last = TRUE)
  return(ret)
}


do_genewise_qc <- function(expression_vector, all_nt_idxs, control_group_complement, curr_df, grna_group_idxs, n_nonzero_trt_thresh, n_nonzero_cntrl_thresh) {
  # 1. compute n nonzero treatment per grna group
  grna_group_posits <- match(x = curr_df$grna_group, table = names(grna_group_idxs))
  n_nonzero_trt_curr <- compute_n_nonzero_trt_vector(expression_vector = expression_vector,
                                                     grna_group_idxs = grna_group_idxs,
                                                     grna_group_posits = grna_group_posits)
  # 2. compute n nonzero control
  if (!control_group_complement) { # nt cell control group
    expression_vector_nt <- expression_vector[all_nt_idxs]
    n_nonzero_cntrl_curr <- sum(expression_vector_nt >= 1)
  } else { # complement set control group
    n_nonzero_total <- sum(expression_vector >= 1)
    n_nonzero_cntrl_curr <- n_nonzero_total - n_nonzero_trt_curr
  }

  # 3. put these vectors into the curr_df
  curr_df$n_nonzero_trt <- n_nonzero_trt_curr
  curr_df$n_nonzero_cntrl <- n_nonzero_cntrl_curr

  # 4. determine the rows of curr_df that pass pairwise qc
  pass_qc <- (n_nonzero_trt_curr >= n_nonzero_trt_thresh) & (n_nonzero_cntrl_curr >= n_nonzero_cntrl_thresh)

  # 5. check if any rows pass pairwise qc
  any_pass_qc <- any(pass_qc)

  # 6. return a list containing pass_qc and if any pass qc
  return(list(curr_df = curr_df, pass_qc = pass_qc, any_pass_qc = any_pass_qc))
}


run_crt_in_memory <- function(response_matrix, grna_assignments,
                              covariate_matrix, response_grna_group_pairs,
                              return_resampling_dist, fit_skew_normal,
                              B1, B2, B3, calibration_check, control_group, n_nonzero_trt_thresh,
                              n_nonzero_cntrl_thresh, return_debugging_metrics, print_progress) {
  if (!control_group_complement) stop("Function not yet implemented for control group NT cells.")
  # 1. setup
  result_list_outer <- vector(mode = "list", length = 2 * length(unique(response_grna_group_pairs$response_id)))
  precomp_list <- vector(mode = "list", length(unique(response_grna_group_pairs$response_id)))
  out_counter <- 1L
  all_nt_idxs <- if (!is.null(grna_assignments$all_nt_idxs)) grna_assignments$all_nt_idxs else NA
  grna_group_idxs <- if (!is.null(grna_assignments$grna_group_idxs)) grna_assignments$grna_group_idxs else NA
  indiv_nt_grna_idxs <- if (!is.null(grna_assignments$indiv_nt_grna_idxs)) grna_assignments$indiv_nt_grna_idxs else NA
  if (calibration_check) {
    low_level_association_funct <- "calibration_complement_crt"
  } else {
    low_level_association_funct <- "discovery_complement_crt"
  }
  n_cells <- nrow(covariate_matrix)

  # 2. perform the gene precomputations; loop over each gene and regress that gene onto the technical factors
  response_ids <- unique(response_grna_group_pairs$response_id)
  for (response_idx in seq_along(response_ids)) {
    if ((response_idx == 1 || response_idx %% 10 == 0) && print_progress) {
      cat(paste0("Running precomputation on gene ", as.character(response_ids[response_idx]), " (", response_idx, " of ", length(response_ids), ")\n"))
    }
    if (response_idx %% 200 == 0) gc() |> invisible()
    response_id <- as.character(response_ids[response_idx])
    # 3. load the expressions of the current response
    expression_vector <- load_csr_row(j = response_matrix@j,
                                      p = response_matrix@p,
                                      x = response_matrix@x,
                                      row_idx = which(rownames(response_matrix) == response_id),
                                      n_cells = n_cells)

    # 4. perform the gene precomputation and add to precomputation list
    response_precomp <- perform_response_precomputation(expressions = expression_vector,
                                                        covariate_matrix = covariate_matrix,
                                                        regression_method = "pois_glm")
    precomp_list[[response_idx]] <- response_precomp

    # 5. perform the gene-wise QC

  }
  names(precomp_list) <- response_ids

  # 5. loop over gRNA groups
  grna_groups <- unique(response_grna_group_pairs$grna_group)
  for (grna_group_idx in seq_along(grna_groups)) {
    if ((grna_group_idx == 1 || grna_group_idx %% 5 == 0) && print_progress) {
      cat(paste0("Analyzing pairs containing gRNA group ", as.character(grna_groups[grna_group_idx]), " (", grna_group_idx, " of ", length(grna_groups), ")\n"))
    }
    if (grna_group_idx %% 200 == 0) gc() |> invisible()
    curr_grna_group <- as.character(grna_groups[grna_group_idx])

    # 6. obtain the set of indices for this grna group (consider making into a function)
    idxs <- if (calibration_check) {
      get_idx_vector_calibration_check(curr_grna_group = curr_grna_group,
                                       indiv_nt_grna_idxs = indiv_nt_grna_idxs,
                                       take_unique = FALSE)
    } else {
      get_idx_vector_discovery_analysis(curr_grna_group = curr_grna_group,
                                        grna_group_idxs = grna_group_idxs)
    }

    trt_idxs <- idxs$trt_idxs
    n_trt <- idxs$n_trt

    # 7. perform the grna precomputation
    fitted_probabilities <- perform_grna_precomputation(trt_idxs = idxs$trt_idxs,
                                                        covariate_matrix = covariate_matrix,
                                                        return_fitted_values = TRUE)

    # 8. obtain the synthetic grna group indices (consider also other idea for sampling the crt indices: fix cell; get probability for that cell; draw single binomial sample s with that probability across B; then, sample WOR s elements from [1, ..., B]. Finally, do a sparse matrix transpose operation. Also consider Gene's idea.)
    synthetic_idxs <- crt_index_sampler_fast(fitted_probabilities = fitted_probabilities,
                                             B = B1 + B2)

    # 9. obtain the grna groups to analyze
    l <- response_grna_group_pairs$grna_group == curr_grna_group
    curr_df <- response_grna_group_pairs[l,]

    ####################################
    # Possibly inside low-level function
    ####################################

    # 10. loop over the genes
    curr_response_ids <- as.character(curr_df$response_id)
    result_list_inner <- vector(mode = "list", length = length(curr_response_ids))
    for (curr_response_idx in seq_along(curr_response_ids)) {
      curr_response_id <- curr_response_ids[curr_response_idx]
      # obtain precomputation for gene
      curr_expression_vector <- load_csr_row(j = response_matrix@j,
                                             p = response_matrix@p,
                                             x = response_matrix@x,
                                             row_idx = which(rownames(response_matrix) == curr_response_id),
                                             n_cells = n_cells)
      response_precomp <- precomp_list[[curr_response_id]]

      # get the precomputation pieces (consider option to not recompute D)
      precomp_pieces <- compute_precomputation_pieces(expression_vector = curr_expression_vector,
                                                      covariate_matrix = covariate_matrix,
                                                      fitted_coefs = response_precomp$fitted_coefs,
                                                      theta = response_precomp$theta,
                                                      full_test_stat = TRUE)

      # compute the result
      result <- run_low_level_test_full_v4(y = curr_expression_vector,
                                           mu = precomp_pieces$mu,
                                           a = precomp_pieces$a,
                                           w = precomp_pieces$w,
                                           D = precomp_pieces$D,
                                           n_trt = n_trt,
                                           use_all_cells = TRUE,
                                           trt_idxs = idxs$trt_idxs,
                                           synthetic_idxs = synthetic_idxs,
                                           B1 = B1, B2 = B2, B3 = 0L,
                                           fit_skew_normal = fit_skew_normal,
                                           return_resampling_dist = return_resampling_dist)
      result_list_inner[[curr_response_idx]] <- result
    }
    # combine the grna group-wise results
    result_list_outer[[out_counter]] <- construct_data_frame_v2(curr_df, result_list_inner,
                                                                return_debugging_metrics, return_resampling_dist)
    out_counter <- out_counter + 1L
  }

  # combine and sort result
  ret <- data.table::rbindlist(result_list_outer, fill = TRUE)
  data.table::setorderv(ret, cols = c("p_value", "response_id"), na.last = TRUE)
  return(ret)
}
