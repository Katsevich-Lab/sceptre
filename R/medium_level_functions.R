run_lowmoi_in_memory <- function(response_matrix, grna_assignments,
                                 covariate_matrix, response_grna_group_pairs,
                                 synthetic_idxs, full_test_stat, return_resampling_dist,
                                 fit_skew_normal, B1, B2, B3, calibration_check,
                                 n_nonzero_trt_thresh, n_nonzero_cntrl_thresh,
                                 return_debugging_metrics, print_progress) {
  # 0. preliminary setup; initialize the args_to_pass, set the low_level_association_funct
  result_list_outer <- vector(mode = "list", length = 2 * length(unique(response_grna_group_pairs$response_id)))
  out_counter <- 1L
  args_to_pass <- list(synthetic_idxs = synthetic_idxs,
                       B1 = B1, B2 = B2, B3 = B3,
                       fit_skew_normal = fit_skew_normal,
                       return_resampling_dist = return_resampling_dist)
  if (calibration_check) {
    args_to_pass$indiv_nt_grna_idxs <- grna_assignments$indiv_nt_grna_idxs
  } else {
    args_to_pass$grna_group_idxs <- grna_assignments$grna_group_idxs
    args_to_pass$covariate_matrix <- covariate_matrix
  }

  low_level_association_funct <- if (!calibration_check & full_test_stat) {
    "lowmoi_full_stat_discovery"
  } else if (!calibration_check & !full_test_stat) {
    "lowmoi_distilled_stat_discovery"
  } else if (calibration_check & full_test_stat) {
    "lowmoi_full_stat_undercover"
  } else {
    "lowmoi_distilled_stat_undercover"
  }

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

    # 3. load the expressions of the current response
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
    response_precomp <- perform_response_precomputation(expressions = expression_vector_nt,
                                                        covariate_matrix = covariate_matrix_nt)

    # 7. obtain precomputation peices for NT cells
    pieces_precomp <- compute_precomputation_pieces(expression_vector_nt,
                                                    covariate_matrix_nt,
                                                    response_precomp$fitted_coefs,
                                                    response_precomp$theta,
                                                    full_test_stat)

    # 8. update the args to pass with grna_groups, expression_vector, response_precomp, pieces_precomp
    args_to_pass$grna_groups <- as.character(curr_df$grna_group)
    args_to_pass$pieces_precomp <- pieces_precomp
    args_to_pass$expression_vector_nt <- expression_vector_nt
    if (!calibration_check) {
      args_to_pass$expression_vector <- expression_vector
      args_to_pass$response_precomp <- response_precomp
    }

    # 9. pass the arguments to the appropriate low-level association testing function
    curr_response_result <- do.call(what = low_level_association_funct, args = args_to_pass)

    # 10. combine the response-wise results into a data table; insert into list
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
