compute_pairwise_qc_information <- function(sceptre_object) {
  # extract basic information from sceptre_object
  cells_in_use <- sceptre_object@cells_in_use
  grna_assignments <- sceptre_object@grna_assignments
  response_matrix <- get_response_matrix(sceptre_object)
  discovery_pairs <- sceptre_object@discovery_pairs
  positive_control_pairs <- sceptre_object@positive_control_pairs
  control_group_complement <- sceptre_object@control_group_complement
  treatment_group_inclusive <- sceptre_object@treatment_group_inclusive
  grna_group_idxs <- grna_assignments$grna_group_idxs

  # 1. construct the entire set of response pairs by combining discovery pairs and pc pairs
  response_grna_group_pairs <- rbind(
    dplyr::mutate(discovery_pairs, discovery = TRUE),
    dplyr::mutate(positive_control_pairs, discovery = FALSE)
  )

  # 2. map each response id and grna group to its index
  dt <- data.table::data.table(
    response_idx = match(
      x = response_grna_group_pairs$response_id,
      table = rownames(response_matrix)
    ),
    response_id = response_grna_group_pairs$response_id,
    grna_idx = match(
      x = response_grna_group_pairs$grna_group,
      table = names(grna_group_idxs)
    ),
    grna_group = response_grna_group_pairs$grna_group,
    discovery = response_grna_group_pairs$discovery
  ) |>
    data.table::setorder(cols = "response_idx")
  to_analyze_response_idxs <- dt$response_idx
  to_analyze_grna_idxs <- dt$grna_idx

  # take cases on odm
  if (methods::is(response_matrix, "odm")) {
    if (sceptre_object@nuclear) {
      out <- ondisc:::compute_n_ok_pairs_ondisc(
        file_name_in = response_matrix@h5_file,
        f_row_ptr = response_matrix@ptr,
        n_genes = nrow(response_matrix),
        n_cells_orig = ncol(response_matrix),
        n_cells_sub = length(cells_in_use),
        grna_group_idxs = grna_group_idxs,
        all_nt_idxs = if (!control_group_complement) grna_assignments$all_nt_idxs else integer(),
        to_analyze_response_idxs = to_analyze_response_idxs,
        to_analyze_grna_idxs = to_analyze_grna_idxs,
        control_group_complement = control_group_complement,
        cells_in_use = cells_in_use,
        unique_response_idxs = sort(unique(to_analyze_response_idxs))
      )
      out$n_nonzero_mat <- matrix()
      out$n_nonzero_tot <- integer()
    } else {
      out <- ondisc:::compute_nt_nonzero_matrix_and_n_ok_pairs_ondisc(
        file_name_in = response_matrix@h5_file,
        f_row_ptr = response_matrix@ptr,
        n_genes = nrow(response_matrix),
        n_cells_orig = ncol(response_matrix),
        n_cells_sub = length(cells_in_use),
        grna_group_idxs = grna_group_idxs,
        indiv_nt_grna_idxs = grna_assignments$indiv_nt_grna_idxs,
        all_nt_idxs = if (!control_group_complement) grna_assignments$all_nt_idxs else integer(),
        to_analyze_response_idxs = to_analyze_response_idxs,
        to_analyze_grna_idxs = to_analyze_grna_idxs,
        control_group_complement = control_group_complement,
        cells_in_use = cells_in_use
      )
    }
  } else {
    out <- compute_nt_nonzero_matrix_and_n_ok_pairs_v3(
      j = response_matrix@j,
      p = response_matrix@p,
      n_cells_orig = ncol(response_matrix),
      n_cells_sub = length(cells_in_use),
      grna_group_idxs = grna_group_idxs,
      indiv_nt_grna_idxs = grna_assignments$indiv_nt_grna_idxs,
      all_nt_idxs = if (!control_group_complement) grna_assignments$all_nt_idxs else integer(),
      to_analyze_response_idxs = to_analyze_response_idxs,
      to_analyze_grna_idxs = to_analyze_grna_idxs,
      control_group_complement = control_group_complement,
      treatment_group_inclusive = treatment_group_inclusive,
      cells_in_use = cells_in_use
    )
  }

  # process results and update sceptre_object
  dt$n_nonzero_trt <- out$n_nonzero_trt
  dt$n_nonzero_cntrl <- out$n_nonzero_cntrl
  sceptre_object@discovery_pairs_with_info <- dt[dt$discovery, c("response_id", "grna_group", "n_nonzero_trt", "n_nonzero_cntrl")]
  sceptre_object@positive_control_pairs_with_info <- dt[!dt$discovery, c("response_id", "grna_group", "n_nonzero_trt", "n_nonzero_cntrl")]
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
    if (nrow(methods::slot(sceptre_object, data_frame_name)) >= 1L) {
      methods::slot(sceptre_object, data_frame_name) <- methods::slot(sceptre_object, data_frame_name) |>
        dplyr::mutate(pass_qc = (n_nonzero_trt >= n_nonzero_trt_thresh & n_nonzero_cntrl >= n_nonzero_cntrl_thresh))
      if (sceptre_object@grna_integration_strategy == "bonferroni") {
        grna_target_data_frame <- sceptre_object@grna_target_data_frame |>
          dplyr::select(grna_id, grna_target)
        n_ok_pairs <- methods::slot(sceptre_object, data_frame_name) |>
          dplyr::select("grna_id" = "grna_group", pass_qc, response_id) |>
          dplyr::left_join(grna_target_data_frame, by = "grna_id") |>
          dplyr::group_by(response_id, grna_target) |>
          dplyr::summarize(any_pass_qc = any(pass_qc), .groups = "drop") |>
          dplyr::pull(any_pass_qc) |>
          sum()
      } else {
        n_ok_pairs <- sum(methods::slot(sceptre_object, data_frame_name)$pass_qc)
      }
      methods::slot(sceptre_object, n_ok_pairs_name) <- n_ok_pairs
    }
  }
  return(sceptre_object)
}
