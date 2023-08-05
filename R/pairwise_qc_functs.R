compute_pairwise_qc_information <- function(sceptre_object) {
  grna_assignments <- sceptre_object@grna_assignments
  response_matrix <- sceptre_object@response_matrix
  discovery_pairs <- sceptre_object@discovery_pairs
  positive_control_pairs <- sceptre_object@positive_control_pairs
  control_group_complement <- sceptre_object@control_group_complement
  compute_effective_sample_sizes <- (nrow(discovery_pairs) >= 1L || nrow(positive_control_pairs) >= 1L)
  if (compute_effective_sample_sizes) {
    grna_group_idxs <- grna_assignments$grna_group_idxs
    response_grna_group_pairs <- rbind(dplyr::mutate(discovery_pairs, discovery = TRUE),
                                       dplyr::mutate(positive_control_pairs, discovery = FALSE))
    dt <- data.table::data.table(response_idx = match(x = response_grna_group_pairs$response_id,
                                                      table = rownames(response_matrix)),
                                 response_id = response_grna_group_pairs$response_id,
                                 grna_idx = match(x = response_grna_group_pairs$grna_group,
                                                  table = names(grna_group_idxs)),
                                 grna_group = response_grna_group_pairs$grna_group,
                                 discovery = response_grna_group_pairs$discovery) |>
      data.table::setorder(cols = "response_idx")
    to_analyze_response_idxs <- dt$response_idx
    to_analyze_grna_idxs <- dt$grna_idx
  } else {
    to_analyze_response_idxs <- to_analyze_grna_idxs <- integer()
    grna_group_idxs <- list()
  }
  out <- compute_nt_nonzero_matrix_and_n_ok_pairs_v2(j = response_matrix@j,
                                                     p = response_matrix@p,
                                                     n_cells = ncol(response_matrix),
                                                     grna_group_idxs = grna_group_idxs,
                                                     indiv_nt_grna_idxs = grna_assignments$indiv_nt_grna_idxs,
                                                     all_nt_idxs = if (!control_group_complement) grna_assignments$all_nt_idxs else integer(),
                                                     to_analyze_response_idxs = to_analyze_response_idxs,
                                                     to_analyze_grna_idxs = to_analyze_grna_idxs,
                                                     compute_n_ok_pairs = compute_effective_sample_sizes,
                                                     control_group_complement = control_group_complement)
  if (compute_effective_sample_sizes) {
    dt$n_nonzero_trt <- out$n_nonzero_trt
    dt$n_nonzero_cntrl <- out$n_nonzero_cntrl
    sceptre_object@discovery_pairs_with_info <- dt[dt$discovery, c("response_id", "grna_group", "n_nonzero_trt", "n_nonzero_cntrl")]
    sceptre_object@positive_control_pairs_with_info <- dt[!dt$discovery, c("response_id", "grna_group", "n_nonzero_trt", "n_nonzero_cntrl")]
  }

  # update the fields of the sceptre_object
  sceptre_object@M_matrix <- out$n_nonzero_mat
  sceptre_object@n_nonzero_tot_vector <- out$n_nonzero_tot
  return(sceptre_object)
}


compute_qc_metrics <- function(sceptre_object) {
  n_nonzero_trt_thresh <- sceptre_object@n_nonzero_trt_thresh
  n_nonzero_cntrl_thresh <- sceptre_object@n_nonzero_cntrl_thresh
  data_frame_names <- c("discovery_pairs_with_info", "positive_control_pairs_with_info")
  n_ok_pair_names <- c("n_ok_discovery_pairs", "n_ok_positive_control_pairs")
  for (i in c(1L, 2L)) {
    data_frame_name <- data_frame_names[i]
    n_ok_pairs_name <- n_ok_pair_names[i]
    if (nrow(slot(sceptre_object, data_frame_name)) >= 1L) {
      slot(sceptre_object, data_frame_name) <- slot(sceptre_object, data_frame_name) |>
        dplyr::mutate(pass_qc = (n_nonzero_trt >= n_nonzero_trt_thresh & n_nonzero_cntrl >= n_nonzero_cntrl_thresh))
      slot(sceptre_object, n_ok_pairs_name) <- sum(slot(sceptre_object, data_frame_name)$pass_qc)
    }
  }
  return(sceptre_object)
}
