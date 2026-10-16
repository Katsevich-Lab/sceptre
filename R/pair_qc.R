# Measurement and policy are independent of the S4 facade. The storage adapter
# retains its fused NT summary scan; policy consumes only aligned count vectors.
#' Collect pair counts without applying a QC decision
#'
#' @inheritParams storage_pair_counts
#' @return The storage adapter's count list, with identical analysis/source index
#'   contracts, request order, optional summary shape and ownership guarantees.
#' @details Keeps the fused NT summary scan in storage; performs no S4 mutation,
#'   threshold policy, nuisance fitting, RNG draws or file writes.
#' @noRd
collect_pair_qc_metrics <- function(store, pair_index, groups, source_cells,
                                    controls, nt_groups = NULL, nuclear = FALSE) {
    storage_pair_counts(store, pair_index, groups, source_cells,
                        controls, nt_groups, nuclear)
}

#' Apply the shared observed-count eligibility rule in pair order
#'
#' @param metrics List of aligned integer `n_nonzero_trt`, `n_nonzero_cntrl`
#'   vectors, as returned by the storage adapter.
#' @param qc_options List with validated nonmissing integer
#'   `n_nonzero_trt_thresh`, `n_nonzero_cntrl_thresh` scalars.
#' @return Logical vector, one per pair: both counts must meet their threshold.
#'   Missing counts preserve the old R `&` semantics (a known failure wins over
#'   an unknown count); unequal count-vector lengths error and empty inputs return
#'   `logical()`. The native calibration sampler shares the same predicate.
#' @details Does not alter measurements, mutate inputs, draw RNG or access files.
#' @noRd
evaluate_pair_eligibility <- function(metrics, qc_options) {
    evaluate_pair_eligibility_cpp(
        metrics$n_nonzero_trt, metrics$n_nonzero_cntrl,
        qc_options$n_nonzero_trt_thresh, qc_options$n_nonzero_cntrl_thresh
    )
}
