#' Add treatment and control counts to discovery and positive-control pairs
#'
#' Prepare cell/group lookups, count requested pairs through the storage
#' helpers,
#' and save the counts in the object's discovery and positive-control tables.
#' Also save the NT-group-by-gene summary matrix and per-gene totals needed for
#' calibration. Keep each table's requested pairs aligned with its returned
#' counts.
#' @noRd
compute_pairwise_qc_information <- function(sceptre_object) {
    # extract basic information from sceptre_object
    cells_in_use <- sceptre_object@preparation@analysis_cells
    grna_assignments <- sceptre_object@preparation@analysis_groups
    response_matrix <- get_response_matrix(sceptre_object)
    discovery_pairs <- sceptre_object@preparation@discovery_pairs
    positive_control_pairs <- sceptre_object@preparation@positive_control_pairs
    control_group_complement <- sceptre_object@spec@control_group_complement
    grna_group_idxs <- grna_assignments$grna_group_idxs

    response_grna_group_pairs <- rbind(
        dplyr::mutate(discovery_pairs, discovery = TRUE),
        dplyr::mutate(positive_control_pairs, discovery = FALSE)
    )

    # Use the same link between original matrix columns and cells kept after QC
    # as the association tests, so counts refer to the same cells.
    index <- prepare_index_context(
        cells_in_use,
        ncol(response_matrix),
        grna_assignments,
        response_grna_group_pairs,
        rownames(response_matrix),
        FALSE,
        control_group_complement,
        sceptre_object@data@low_moi
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
        nt_groups = nt_groups
    )

    # process results and update sceptre_object
    dt$n_nonzero_trt <- out$n_nonzero_trt
    dt$n_nonzero_cntrl <- out$n_nonzero_cntrl
    sceptre_object@preparation@discovery_pairs_with_info <- dt[
        dt$discovery,
        c("response_id", "grna_group", "n_nonzero_trt", "n_nonzero_cntrl")
    ]
    sceptre_object@preparation@positive_control_pairs_with_info <- dt[
        !dt$discovery,
        c("response_id", "grna_group", "n_nonzero_trt", "n_nonzero_cntrl")
    ]
    sceptre_object@preparation@nt_nonzero_matrix <- out$n_nonzero_mat
    sceptre_object@preparation@control_nonzero_totals <- out$n_nonzero_tot
    return(sceptre_object)
}


#' Mark QC-passing pairs and count how many hypotheses remain
#' @param sceptre_object Object with treatment/control counts, QC thresholds and
#'   guide integration strategy already set.
#' @return Object with pass_qc in each nonempty pair table and updated
#'   passing-pair
#'   counts. For Bonferroni, count distinct gene-target pairs with any passing
#'   guide.
#'   Otherwise count passing rows, including repeated requests. Leave empty
#'   tables
#'   and their existing count fields unchanged.
#' @noRd
compute_qc_metrics <- function(sceptre_object) {
    n_nonzero_trt_thresh <- sceptre_object@spec@n_nonzero_trt_thresh
    n_nonzero_cntrl_thresh <- sceptre_object@spec@n_nonzero_cntrl_thresh
    data_frame_names <- c(
        "discovery_pairs_with_info",
        "positive_control_pairs_with_info"
    )
    for (i in c(1L, 2L)) {
        data_frame_name <- data_frame_names[i]
        if (
            nrow(methods::slot(sceptre_object@preparation, data_frame_name)) >=
                1L
        ) {
            methods::slot(
                sceptre_object@preparation,
                data_frame_name
            ) <- methods::slot(
                sceptre_object@preparation,
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
        }
    }
    return(sceptre_object)
}
