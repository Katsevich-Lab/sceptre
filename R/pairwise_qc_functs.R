compute_pairwise_qc_information <- function(sceptre_object) {
    # extract basic information from sceptre_object
    cells_in_use <- sceptre_object@cells_in_use
    grna_assignments <- sceptre_object@grna_assignments
    response_matrix <- get_response_matrix(sceptre_object)
    discovery_pairs <- sceptre_object@discovery_pairs
    positive_control_pairs <- sceptre_object@positive_control_pairs
    control_group_complement <- sceptre_object@control_group_complement
    grna_group_idxs <- grna_assignments$grna_group_idxs

    # 1. construct the entire set of response pairs by combining discovery pairs
    # and pc pairs
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
    controls <- if (control_group_complement) {
        "complement"
    } else {
        grna_assignments$all_nt_idxs
    }
    nt_groups <- grna_assignments$indiv_nt_grna_idxs
    if (!control_group_complement) {
        # Legacy objects store NT memberships relative to the ordered NT view.
        # Normalize that persisted representation before crossing the boundary.
        nt_groups <- lapply(nt_groups, function(group) controls[group])
    }
    out <- storage_pair_counts(
        store = response_matrix,
        pair_index = data.frame(response_idx = dt$response_idx, group_idx = dt$grna_idx),
        groups = grna_group_idxs, source_cells = cells_in_use,
        controls = controls, nt_groups = nt_groups,
        nuclear = sceptre_object@nuclear
    )

    # process results and update sceptre_object
    dt$n_nonzero_trt <- out$n_nonzero_trt
    dt$n_nonzero_cntrl <- out$n_nonzero_cntrl
    sceptre_object@discovery_pairs_with_info <- dt[
        dt$discovery,
        c("response_id", "grna_group", "n_nonzero_trt", "n_nonzero_cntrl")
    ]
    sceptre_object@positive_control_pairs_with_info <- dt[
        !dt$discovery,
        c("response_id", "grna_group", "n_nonzero_trt", "n_nonzero_cntrl")
    ]
    sceptre_object@M_matrix <- out$n_nonzero_mat
    sceptre_object@n_nonzero_tot_vector <- out$n_nonzero_tot
    return(sceptre_object)
}


compute_qc_metrics <- function(sceptre_object) {
    n_nonzero_trt_thresh <- sceptre_object@n_nonzero_trt_thresh
    n_nonzero_cntrl_thresh <- sceptre_object@n_nonzero_cntrl_thresh
    data_frame_names <- c(
        "discovery_pairs_with_info",
        "positive_control_pairs_with_info"
    )
    n_ok_pair_names <- c("n_ok_discovery_pairs", "n_ok_positive_control_pairs")
    for (i in c(1L, 2L)) {
        data_frame_name <- data_frame_names[i]
        n_ok_pairs_name <- n_ok_pair_names[i]
        if (nrow(methods::slot(sceptre_object, data_frame_name)) >= 1L) {
            methods::slot(sceptre_object, data_frame_name) <- methods::slot(
                sceptre_object,
                data_frame_name
            ) |>
                dplyr::mutate(
                    pass_qc = (n_nonzero_trt >= n_nonzero_trt_thresh &
                        n_nonzero_cntrl >= n_nonzero_cntrl_thresh)
                )
            if (sceptre_object@grna_integration_strategy == "bonferroni") {
                grna_target_data_frame <- sceptre_object@grna_target_data_frame |>
                    dplyr::select(grna_id, grna_target)
                n_ok_pairs <- methods::slot(sceptre_object, data_frame_name) |>
                    dplyr::select(
                        "grna_id" = "grna_group",
                        pass_qc,
                        response_id
                    ) |>
                    dplyr::left_join(grna_target_data_frame, by = "grna_id") |>
                    dplyr::group_by(response_id, grna_target) |>
                    dplyr::summarize(
                        any_pass_qc = any(pass_qc),
                        .groups = "drop"
                    ) |>
                    dplyr::pull(any_pass_qc) |>
                    sum()
            } else {
                n_ok_pairs <- sum(
                    methods::slot(sceptre_object, data_frame_name)$pass_qc
                )
            }
            methods::slot(sceptre_object, n_ok_pairs_name) <- n_ok_pairs
        }
    }
    return(sceptre_object)
}
