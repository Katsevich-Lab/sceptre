#' Choose negative-control gene-guide pairs that pass the QC thresholds
#' @param sceptre_object Post-QC object containing nontargeting (NT) guide cells,
#'   `M_matrix` (NT-group counts by gene) and `n_nonzero_tot_vector` (totals by gene).
#' @param n_calibration_pairs Requested number of distinct calibration pairs.
#' @param calibration_group_size Number of NT guides to combine per group; must
#'   be positive and no larger than the number of available NT guides.
#' @return Table with gene names, joined NT guide names (`grna_group`), treatment
#'   and control counts, and `pass_qc = TRUE`, in the C++ sampler's order. If too
#'   few pairs qualify, return fewer with a message; if none qualify, error.
#'   C++ sampling uses fixed seeds without advancing R's random-number generator.
#' @noRd
construct_negative_control_pairs <- function(
    sceptre_object,
    n_calibration_pairs,
    calibration_group_size
) {
    message("Constructing negative control pairs.")
    grna_assignments <- sceptre_object@preparation@grna_assignments
    response_matrix <- get_response_matrix(sceptre_object)
    n_nonzero_trt_thresh <- sceptre_object@spec@n_nonzero_trt_thresh
    n_nonzero_cntrl_thresh <- sceptre_object@spec@n_nonzero_cntrl_thresh
    n_nonzero_m <- sceptre_object@preparation@M_matrix
    n_nonzero_tot <- sceptre_object@preparation@n_nonzero_tot_vector
    cells_in_use <- sceptre_object@preparation@cells_in_use
    discovery_pairs_with_info <- sceptre_object@preparation@discovery_pairs_with_info

    N_POSSIBLE_GROUPS_THRESHOLD <- 100L
    nt_grna_names <- names(grna_assignments[["indiv_nt_grna_idxs"]])
    n_nt_grnas <- length(nt_grna_names)
    n_genes <- nrow(response_matrix)

    n_possible_groups <- choose(n_nt_grnas, calibration_group_size)
    if (n_possible_groups <= N_POSSIBLE_GROUPS_THRESHOLD) {
        # iterate over the entire set of combinations
        n_possible_groups <- as.integer(n_possible_groups)
        possible_groups_m <- iterate_over_combinations(
            n_nt_grnas,
            calibration_group_size,
            n_possible_groups
        )
    } else {
        # sample from the set of combinations
        p_hat <- if (sceptre_object@preparation@n_discovery_pairs == 0L) {
            0.1
        } else {
            mean(discovery_pairs_with_info$pass_qc)
        }
        if (!is.finite(p_hat) || p_hat <= 0) {
            # Discovery targets can all fail QC even when NT groups qualify.
            # Use the necessary total-count condition as a proposal rate; the
            # native pair sampler still checks actual treatment/control counts.
            p_hat <- mean(n_nonzero_tot >=
                n_nonzero_trt_thresh + n_nonzero_cntrl_thresh)
            if (!is.finite(p_hat) || p_hat <= 0) {
                stop("No negative control pair passes pairwise QC.")
            }
        }
        possible_groups_m <- sample_calibration_groups(
            calibration_group_size,
            n_calibration_pairs,
            n_possible_groups,
            n_nt_grnas,
            as.numeric(n_genes),
            N_POSSIBLE_GROUPS_THRESHOLD,
            p_hat
        )
    }

    calculate_ess_using_m_matrix <- sceptre_object@data@low_moi ||
        calibration_group_size == 1L
    if (calculate_ess_using_m_matrix) {
        samp <- sample_calibration_pairs(
            n_nonzero_m = n_nonzero_m,
            n_nonzero_tot = n_nonzero_tot,
            possible_groups_m = possible_groups_m,
            n_genes = n_genes,
            n_calibration_pairs = n_calibration_pairs,
            n_nonzero_trt_thresh = n_nonzero_trt_thresh,
            n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh
        )
    } else {
        nt_grna_group_idxs <- lapply(
            X = seq(1L, nrow(possible_groups_m)),
            function(i) {
                nt_grna_idxs <- possible_groups_m[i, ] + 1L
                grna_assignments$indiv_nt_grna_idxs[nt_grna_idxs] |>
                    unlist() |>
                    unique()
            }
        )
        n_nonzero_m <- storage_group_counts(
            response_matrix,
            seq_len(n_genes),
            nt_grna_group_idxs,
            cells_in_use
        )
        dummy_possibly_groups_m <- matrix(
            seq_along(nt_grna_group_idxs) - 1L,
            ncol = 1L
        )
        samp <- sample_calibration_pairs(
            n_nonzero_m = n_nonzero_m,
            n_nonzero_tot = n_nonzero_tot,
            possible_groups_m = dummy_possibly_groups_m,
            n_genes = n_genes,
            n_calibration_pairs = n_calibration_pairs,
            n_nonzero_trt_thresh = n_nonzero_trt_thresh,
            n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh
        )
    }

    if (length(samp$response_idxs) == 0L) {
        stop("No negative control pair passes pairwise QC.")
    }
    response_ids <- rownames(response_matrix)
    increment_group_indices_in_place(possible_groups_m)
    undercover_groups <- get_undercover_group_names(
        possible_groups_m,
        nt_grna_names
    )
    df <- data.frame(
        response_id = response_ids[samp$response_idxs],
        grna_group = undercover_groups[samp$grna_group_idxs],
        n_nonzero_trt = samp$n_nonzero_trt_v,
        n_nonzero_cntrl = samp$n_nonzero_cntrl_v,
        pass_qc = TRUE
    )
    message(crayon::green(" \u2713"))

    if (nrow(df) < n_calibration_pairs) {
        note_text <- paste0(
            "Note: Unable to generate the number of negative control pairs (",
            n_calibration_pairs,
            ") requested. Generating as many negative control pairs (",
            nrow(df),
            ") as possible."
        )
        message(crayon::red(note_text))
    }
    return(df)
}


#' Choose how many NT guides to combine in a calibration group by default
#' @param grna_target_data_frame Guide metadata with `grna_group` already set.
#'   The caller needs at least one targeting group to obtain a finite median.
#' @return Rounded median number of guides per targeting group, limited to the
#'   number of NT metadata rows. This helper adds no minimum or fallback for
#'   missing targeting groups.
#' @noRd
compute_calibration_group_size <- function(grna_target_data_frame) {
    median_group_size <- grna_target_data_frame |>
        dplyr::filter(grna_group != "non-targeting") |>
        dplyr::pull(grna_group) |>
        table() |>
        stats::median() |>
        round()
    n_ntc <- grna_target_data_frame |>
        dplyr::filter(grna_group == "non-targeting") |>
        nrow()
    calibration_group_size <- as.integer(min(median_group_size, n_ntc))
    return(calibration_group_size)
}


#' Join the NT guide names that make up each calibration group
#' @param possible_groups_m Integer matrix with one guide combination per row.
#'   Entries are positions in `nt_grna_names`, starting at 1 after the C++ output
#'   has been converted to R numbering.
#' @param nt_grna_names NT guide names in the same order as the assignment list.
#' @return One name per matrix row, for example `"NT1&NT2"`, keeping row order
#'   and guide order within each combination.
#' @noRd
get_undercover_group_names <- function(possible_groups_m, nt_grna_names) {
    apply(X = possible_groups_m, MARGIN = 1, FUN = function(r) {
        paste0(nt_grna_names[r], collapse = "&")
    })
}
