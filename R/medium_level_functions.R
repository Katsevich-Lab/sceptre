run_lowmoi_in_memory <- function(response_matrix, grna_assignments,
                                 covariate_matrix, response_grna_group_pairs,
                                 synthetic_idxs, return_resampling_dist, fit_skew_normal,
                                 B1, B2, B3, calibration_check, discovery_test_stat,
                                 n_nonzero_trt_thresh, n_nonzero_cntrl_thresh,
                                 return_debugging_metrics, regression_method, print_progress) {
  # 0. preliminary setup; initialize the args_to_pass, set the low_level_association_funct
  result_list_outer <- vector(mode = "list", length = 2 * length(unique(response_grna_group_pairs$response_id)))
  out_counter <- 1L
  args_to_pass <- list(synthetic_idxs = synthetic_idxs, B1 = B1, B2 = B2, B3 = B3,
                       fit_skew_normal = fit_skew_normal,
                       return_resampling_dist = return_resampling_dist,
                       grna_group_idxs = grna_assignments$grna_group_idxs,
                       covariate_matrix = covariate_matrix,
                       all_nt_idxs = grna_assignments$all_nt_idxs,
                       regression_method = regression_method,
                       indiv_nt_grna_idxs = if (calibration_check) grna_assignments$indiv_nt_grna_idxs else NA)
  SE_THRESH <- 15.0
  if (calibration_check) {
    low_level_association_funct <- "lowmoi_undercover_stat"
  } else {
    low_level_association_funct <- if (discovery_test_stat == "exact") "lowmoi_exact_stat_discovery" else "lowmoi_approximate_stat_discovery"
  }
  run_outer_regression <- calibration_check || (discovery_test_stat == "approximate")

  # 1. obtain the subset of the covariate matrix corresponding to the NT cells and n_cells
  covariate_matrix_nt <- covariate_matrix[grna_assignments$all_nt_idxs,]
  n_cells <- ncol(response_matrix)

  # 2. loop over the response IDs
  response_ids <- unique(response_grna_group_pairs$response_id)
  for (response_idx in seq_along(response_ids)) {
    if ((response_idx == 1 || response_idx %% 5 == 0) && print_progress) {
      cat(paste0("Analyzing pairs containing response ", as.character(response_ids[response_idx]), " (", response_idx, " of ", length(response_ids), ")\n"))
    }
    if (response_idx %% 200 == 0) gc() |> invisible()
    response_id <- as.character(response_ids[response_idx])

    # 3. load the expressions of the current response; also get the nt expression vector
    expression_vector <- load_csr_row(j = response_matrix@j,
                                      p = response_matrix@p,
                                      x = response_matrix@x,
                                      row_idx = which(rownames(response_matrix) == response_id),
                                      n_cells = n_cells)
    expression_vector_nt <- expression_vector[grna_assignments$all_nt_idxs]

    # 4. obtain the gRNA groups to analyze
    l <- response_grna_group_pairs$response_id == response_id
    curr_df <- response_grna_group_pairs[l,]

    # 5. if running a discovery analysis, do QC
    if (!calibration_check) {
      n_nonzero_cntrl_curr <- sum(expression_vector_nt >= 1)
      grna_group_posits <- match(x = curr_df$grna_group, table = names(grna_assignments$grna_group_idxs))
      n_nonzero_trt_curr <- compute_n_nonzero_trt_vector(expression_vector = expression_vector,
                                                         grna_group_idxs = grna_assignments$grna_group_idxs,
                                                         grna_group_posits = grna_group_posits)
      curr_df$n_nonzero_trt <- n_nonzero_trt_curr
      curr_df$n_nonzero_cntrl <- n_nonzero_cntrl_curr
      pass_qc <- n_nonzero_trt_curr >= n_nonzero_trt_thresh

      # i. if n_nonzero_cntrl_curr is less than n_nonzero_cntrl_thresh, jump to next iteration
      if (n_nonzero_cntrl_curr < n_nonzero_cntrl_thresh || !any(pass_qc)) {
        result_list_outer[[out_counter]] <- curr_df
        out_counter <- out_counter + 1L
        next
      }

      # ii. remove any rows that have not passed qc; keep the rows that have passed qc
      if (!all(pass_qc)) {
        result_list_outer[[out_counter]] <- curr_df[!pass_qc,]
        out_counter <- out_counter + 1L
      }
      curr_df <- curr_df[pass_qc,]
    }


    # 6. perform the expression on technical factor regression
    if (run_outer_regression) {
      response_precomp <- perform_response_precomputation(expressions = expression_vector_nt,
                                                          covariate_matrix = covariate_matrix_nt,
                                                          regression_method = regression_method)

      # 7. obtain precomputation peices for NT cells
      pieces_precomp <- compute_precomputation_pieces(expression_vector_nt,
                                                      covariate_matrix_nt,
                                                      response_precomp$fitted_coefs,
                                                      response_precomp$theta,
                                                      full_test_stat = TRUE)
    } else {
      response_precomp <- pieces_precomp <- NA
    }

    args_to_pass$grna_groups <- as.character(curr_df$grna_group)
    args_to_pass$pieces_precomp <- pieces_precomp
    args_to_pass$expression_vector_nt <- expression_vector_nt
    args_to_pass$expression_vector <- expression_vector
    args_to_pass$response_precomp <- response_precomp

    # 9. pass the arguments to the appropriate low-level association testing function
    curr_response_result <- do.call(what = low_level_association_funct, args = args_to_pass)

    # 10. combine the response-wise results into a data table; insert into list
    result_list_outer[[out_counter]] <- construct_data_frame_v2(curr_df, curr_response_result,
                                                                return_debugging_metrics, return_resampling_dist,
                                                                response_precomp$precomp_str, low_level_association_funct)
    out_counter <- out_counter + 1L
  }

  # combine and sort result
  ret <- data.table::rbindlist(result_list_outer, fill = TRUE)
  data.table::setorderv(ret, cols = c("p_value", "response_id"), na.last = TRUE)
  return(ret)
}
