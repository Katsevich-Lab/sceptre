#' Measure pair sparsity while retaining the saved QC slot layout
#'
#' @param sceptre_object Cellwise-QC object with source-cell `cells_in_use`,
#'   post-QC group memberships and discovery/positive-control pair tables.
#' @return Updated object with pair tables sorted by source response row,
#'   keeping
#'   repeated requests. Counts refer to stored sparse positions, including
#'   explicit
#'   zeros. `M_matrix` has NT-guide rows and all-response columns; totals use
#'   the
#'   configured control view. The nuclear backed path leaves NT summaries empty.
#' @noRd
compute_pairwise_qc_information <- function(sceptre_object) {
    # extract basic information from sceptre_object
    cells_in_use <- sceptre_object@cells_in_use
    grna_assignments <- sceptre_object@grna_assignments
    response_matrix <- get_response_matrix(sceptre_object)
    discovery_pairs <- sceptre_object@discovery_pairs
    positive_control_pairs <- sceptre_object@positive_control_pairs
    control_group_complement <- sceptre_object@control_group_complement
    grna_group_idxs <- grna_assignments$grna_group_idxs

    response_grna_group_pairs <- rbind(
        dplyr::mutate(discovery_pairs, discovery = TRUE),
        dplyr::mutate(positive_control_pairs, discovery = FALSE)
    )

    # Resolve the same source/analysis coordinates used by association workers.
    index <- prepare_index_context(
        cells_in_use,
        ncol(response_matrix),
        grna_assignments,
        response_grna_group_pairs,
        rownames(response_matrix),
        FALSE,
        control_group_complement,
        sceptre_object@low_moi
    )
    dt <- data.table::data.table(
        response_idx = index$pair_rows$response_idx,
        response_id = response_grna_group_pairs$response_id,
        grna_idx = index$pair_rows$group_idx,
        grna_group = response_grna_group_pairs$grna_group,
        discovery = response_grna_group_pairs$discovery
    ) |>
        data.table::setorder(cols = "response_idx")
    controls <- if (control_group_complement) {
        "complement"
    } else {
        index$nt_to_analysis
    }
    nt_groups <- index$nt_groups_analysis
    out <- collect_pair_qc_metrics(
        store = response_matrix,
        pair_index = data.frame(
            response_idx = dt$response_idx,
            group_idx = dt$grna_idx
        ),
        groups = grna_group_idxs,
        source_cells = cells_in_use,
        controls = controls,
        nt_groups = nt_groups,
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


#' Apply pair eligibility thresholds and count retained hypotheses
#'
#' @param sceptre_object Object whose QC tables contain aligned
#'   treatment/control
#'   counts and whose threshold slots and integration strategy have been set.
#' @return Updated object with `pass_qc` in existing row order and retained-pair
#'   counts. Bonferroni counts distinct response/target groups with any passing
#'   guide; other strategies count passing rows, including repeated requests.
#'   Empty tables and their existing count slots are left unchanged.
#' @noRd
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
                    pass_qc = evaluate_pair_eligibility(
                        list(
                            n_nonzero_trt = n_nonzero_trt,
                            n_nonzero_cntrl = n_nonzero_cntrl
                        ),
                        list(
                            n_nonzero_trt_thresh = n_nonzero_trt_thresh,
                            n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh
                        )
                    )
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
