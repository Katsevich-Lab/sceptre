run_perm_test_in_memory <- function(response_matrix, grna_assignments,
                                    covariate_matrix, response_grna_group_pairs,
                                    synthetic_idxs, return_resampling_dist, fit_skew_normal,
                                    B1, B2, B3, calibration_check, control_group, n_nonzero_trt_thresh,
                                    n_nonzero_cntrl_thresh, return_debugging_metrics, print_progress) {
  # 1. setup
  result_list_outer <- vector(mode = "list", length = 2 * length(unique(response_grna_group_pairs$response_id)))
  out_counter <- 1L
  all_nt_idxs <- if (!is.null(grna_assignments$all_nt_idxs)) grna_assignments$all_nt_idxs else NA
  grna_group_idxs <- if (!is.null(grna_assignments$grna_group_idxs)) grna_assignments$grna_group_idxs else NA
  indiv_nt_grna_idxs <- if (!is.null(grna_assignments$indiv_nt_grna_idxs)) grna_assignments$indiv_nt_grna_idxs else NA
  args_to_pass <- list(synthetic_idxs = synthetic_idxs,
                       B1 = B1, B2 = B2, B3 = B3,
                       fit_skew_normal = fit_skew_normal,
                       return_resampling_dist = return_resampling_dist,
                       covariate_matrix = covariate_matrix,
                       all_nt_idxs = all_nt_idxs,
                       grna_group_idxs = grna_group_idxs,
                       indiv_nt_grna_idxs = indiv_nt_grna_idxs)
  if (calibration_check && !control_group_complement) {
    low_level_association_funct <- "calibration_ntcells_perm_test"
  } else if (calibration_check && control_group_complement) {
    low_level_association_funct <- "calibration_complement_perm_test"
  } else if (!calibration_check && !control_group_complement) {
    low_level_association_funct <- "discovery_ntcells_perm_test"
  } else if (!calibration_check && control_group_complement) {
    low_level_association_funct <- "discovery_complement_perm_test"
  }
  run_outer_regression <- low_level_association_funct != "discovery_ntcells_perm_test"
  n_cells <- ncol(response_matrix)
  if (low_level_association_funct == "calibration_ntcells_perm_test") {
    covariate_matrix <- covariate_matrix[all_nt_idxs,]
  }

  # 2. loop over the response IDs
  response_ids <- unique(response_grna_group_pairs$response_id)
  for (response_idx in seq_along(response_ids)) {
    if ((response_idx == 1 || response_idx %% 5 == 0) && print_progress) {
      cat(paste0("Analyzing pairs containing response ", as.character(response_ids[response_idx]), " (", response_idx, " of ", length(response_ids), ")\n"))
    }
    if (response_idx %% 200 == 0) gc() |> invisible()
    response_id <- as.character(response_ids[response_idx])

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
      gene_wise_qc_result <- do_genewise_qc(expression_vector, all_nt_idxs, control_group_complement, curr_df, grna_group_idxs)
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
    } else {
      response_precomp <- pieces_precomp <- NA
    }

    # 7. update the args_to_pass
    args_to_pass$grna_groups <- as.character(curr_df$grna_group)
    args_to_pass$expression_vector <- expression_vector
    args_to_pass$pieces_precomp <- pieces_precomp
    args_to_pass$response_precomp <- response_precomp

    # 8. pass the arguments to the appropriate low-level assocation testing function
    # curr_response_result <- do.call(what = low_level_association_funct, args = args_to_pass)

    # 9. combine the response-wise results into a data table; insert into list
    result_list_outer[[out_counter]] <- construct_data_frame_v2(curr_df, curr_response_result,
                                                                return_debugging_metrics, return_resampling_dist,
                                                                response_precomp$precomp_str)
    out_counter <- out_counter + 1L
  }

  # combine and sort result
  ret <- data.table::rbindlist(result_list_outer, fill = TRUE)
  data.table::setorderv(ret, cols = c("p_value", "response_id"), na.last = TRUE)
  return(ret)
}


do_genewise_qc <- function(expression_vector, all_nt_idxs, control_group_complement, curr_df, grna_group_idxs) {
  # 1. compute n nonzero treatment per grna group
  grna_group_posits <- match(x = curr_df$grna_group, table = names(grna_assignments$grna_group_idxs))
  n_nonzero_trt_curr <- compute_n_nonzero_trt_vector(expression_vector = expression_vector,
                                                     grna_group_idxs = grna_assignments$grna_group_idxs,
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
