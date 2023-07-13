# helper function 1
get_id_from_idx <- function(response_idx, print_progress, response_ids, print_multiple = 5L, gc_multiple = 200L, feature = "response") {
  if ((response_idx == 1 || response_idx %% print_multiple == 0) && print_progress) {
    cat(paste0("Analyzing pairs containing ", feature, " ", as.character(response_ids[response_idx]),
               " (", response_idx, " of ", length(response_ids), ")\n"))
  }
  if (response_idx %% gc_multiple == 0) gc() |> invisible()
  response_id <- as.character(response_ids[response_idx])
  return(response_id)
}

# helper function 2: gene-wise QC
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

# core function 1: run permutation test in memory
run_perm_test_in_memory <- function(response_matrix, grna_assignments, covariate_matrix, response_grna_group_pairs,
                                    synthetic_idxs, output_amount, fit_parametric_curve, B1, B2, B3, calibration_check,
                                    control_group_complement, n_nonzero_trt_thresh, n_nonzero_cntrl_thresh,
                                    side_code, low_moi, print_progress) {
  # 0. define several variables
  gene_pass_qc_list <- gene_fail_qc_list <- vector(mode = "list", length = length(unique(response_grna_group_pairs$response_id)))
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
    if (subset_to_nt_cells) expression_vector <- expression_vector[all_nt_idxs]

    # 4. obtain the gRNA groups to analyze
    l <- response_grna_group_pairs$response_id == response_id
    curr_df <- response_grna_group_pairs[l,]

    if (!calibration_check) {
      gene_wise_qc_result <- do_genewise_qc(expression_vector, all_nt_idxs, control_group_complement,
                                            curr_df, grna_group_idxs, n_nonzero_trt_thresh, n_nonzero_cntrl_thresh)
      gene_fail_qc_list[[response_idx]] <- gene_wise_qc_result$curr_df[!gene_wise_qc_result$pass_qc,]
      if (!gene_wise_qc_result$any_pass_qc) next
      curr_df <- gene_wise_qc_result$curr_df[gene_wise_qc_result$pass_qc,]
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
      curr_response_result <- perm_test_glm_factored_out(synthetic_idxs, B1, B2, B3, fit_parametric_curve,
                                                         output_amount, grna_groups, expression_vector,
                                                         pieces_precomp, get_idx_f, side_code)
    } else {
      curr_response_result <- discovery_ntcells_perm_test(synthetic_idxs, B1, B2, B3, fit_parametric_curve,
                                                          output_amount, covariate_matrix, all_nt_idxs,
                                                          grna_group_idxs, grna_groups, expression_vector, side_code)
    }

    # 9. combine the response-wise results into a data table; insert into list
    gene_pass_qc_list[[response_idx]] <- construct_data_frame_v2(curr_df, curr_response_result, output_amount)
  }


  # combine and sort result
  ret <- data.table::rbindlist(gene_pass_qc_list, fill = TRUE)
  if (!calibration_check) {
    ret_fail_qc <- data.table::rbindlist(gene_fail_qc_list)
    ret <- data.table::rbindlist(list(ret, ret_fail_qc), fill = TRUE)
  }
  data.table::setorderv(ret, cols = c("p_value", "response_id"), na.last = TRUE)
  return(ret)
}


run_crt_in_memory_v2 <- function(response_matrix, grna_assignments, covariate_matrix, response_grna_group_pairs,
                                 output_amount, fit_parametric_curve, B1, B2, B3, calibration_check, control_group_complement,
                                 n_nonzero_trt_thresh, n_nonzero_cntrl_thresh, side_code, low_moi, print_progress) {
  # 0. define several variables
  gene_ess_list <- vector(mode = "list", length = length(unique(response_grna_group_pairs$response_id)))
  subset_to_nt_cells <- calibration_check && !control_group_complement
  run_outer_regression <- calibration_check || control_group_complement
  all_nt_idxs <- if (!is.null(grna_assignments$all_nt_idxs)) grna_assignments$all_nt_idxs else NA
  grna_group_idxs <- if (!is.null(grna_assignments$grna_group_idxs)) grna_assignments$grna_group_idxs else NA
  indiv_nt_grna_idxs <- if (!is.null(grna_assignments$indiv_nt_grna_idxs)) grna_assignments$indiv_nt_grna_idxs else NA
  n_cells <- nrow(covariate_matrix)
  get_idx_f <- get_idx_vector_factory(calibration_check, indiv_nt_grna_idxs, grna_group_idxs, low_moi)
  if (run_outer_regression) gene_precomp_list <- vector(mode = "list", length = length(unique(response_grna_group_pairs$response_id)))

  # 1. subset covariate matrix to nt cells (if applicable)
  if (subset_to_nt_cells) covariate_matrix <- covariate_matrix[all_nt_idxs,]

  # 2. loop over genes, performing the precomputation
  response_ids <- as.character(unique(response_grna_group_pairs$response_id))
  for (response_idx in seq_along(response_ids)) {
    response_id <- as.character(response_ids[response_idx])
    # 3. load the expressions of the current response
    expression_vector <- load_csr_row(j = response_matrix@j,
                                      p = response_matrix@p,
                                      x = response_matrix@x,
                                      row_idx = which(rownames(response_matrix) == response_id),
                                      n_cells = n_cells)
    if (subset_to_nt_cells) expression_vector <- expression_vector[all_nt_idxs]

    # 4. obtain the gRNA groups to analyze
    l <- response_grna_group_pairs$response_id == response_id
    curr_df <- response_grna_group_pairs[l,]

    # 5. run QC if conducting a discovery analysis
    if (!calibration_check) {
      gene_wise_qc_result <- do_genewise_qc(expression_vector, all_nt_idxs, control_group_complement,
                                            curr_df, grna_group_idxs, n_nonzero_trt_thresh, n_nonzero_cntrl_thresh)
      gene_wise_qc_result$curr_df$pass_qc <- gene_wise_qc_result$pass_qc
      gene_ess_list[[response_idx]] <- gene_wise_qc_result$curr_df
      if (!gene_wise_qc_result$any_pass_qc) next
    }

    # 6. perform the gene precomputation and add to precomputation list
    if (run_outer_regression) {
      if ((response_idx == 1 || response_idx %% 5L == 0) && print_progress) {
      cat(paste0("Running precomputation on response ", as.character(response_ids[response_idx]),
                 " (", response_idx, " of ", length(response_ids), ")\n"))
      }
      response_precomp <- perform_response_precomputation(expressions = expression_vector,
                                                          covariate_matrix = covariate_matrix,
                                                          regression_method = "pois_glm")
      gene_precomp_list[[response_idx]] <- response_precomp
    }
  }
  if (run_outer_regression) names(gene_precomp_list) <- response_ids

  # 7. obtain the response-grna group pairs to analyze (if running a discovery analysis)
  if (!calibration_check) {
    response_grna_group_pairs_with_ess <- data.table::rbindlist(gene_ess_list)
    response_grna_group_pairs <- response_grna_group_pairs_with_ess[response_grna_group_pairs_with_ess$pass_qc,]
  }

  # 8. loop over the grna groups
  grna_groups <- unique(response_grna_group_pairs$grna_group)
  result_list_outer <- vector(mode = "list", length = length(grna_groups))
  for (grna_group_idx in seq_along(grna_groups)) {
    curr_grna_group <- get_id_from_idx(grna_group_idx, print_progress, grna_groups, feature = "gRNA group")

    # 11. obtain the genes to analyze
    l <- response_grna_group_pairs$grna_group == curr_grna_group
    curr_df <- response_grna_group_pairs[l,]
    response_ids <- as.character(curr_df$response_id)

    # 12. call the low-level analysis function
    if (run_outer_regression) {
      curr_response_result <- crt_glm_factored_out(B1, B2, fit_parametric_curve, output_amount,
                                                   response_ids, gene_precomp_list, covariate_matrix,
                                                   get_idx_f, curr_grna_group, subset_to_nt_cells, all_nt_idxs,
                                                   n_cells, response_matrix, side_code)
    } else {
      curr_response_result <- discovery_ntcells_crt(B1, B2, fit_parametric_curve, output_amount, get_idx_f,
                                                    response_ids, covariate_matrix, curr_grna_group, all_nt_idxs,
                                                    n_cells, response_matrix, side_code)
    }
    result_list_outer[[grna_group_idx]] <- construct_data_frame_v2(curr_df, curr_response_result, output_amount)
  }
  # combine and sort result
  ret <- data.table::rbindlist(result_list_outer, fill = TRUE)
  if (!calibration_check) {
    ret <- data.table::rbindlist(list(ret,
                                      response_grna_group_pairs_with_ess[
                                        !response_grna_group_pairs_with_ess$pass_qc,]), fill = TRUE)
    ret$pass_qc <- NULL
  }
  data.table::setorderv(ret, cols = c("p_value", "response_id"), na.last = TRUE)
  return(ret)
}
