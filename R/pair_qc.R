# Separate counting from the QC decision. First ask storage for treatment and
# control counts in pair order; then compare those counts with the thresholds.
# Neither helper needs the full sceptre object. Storage can collect the NT
# summaries during the same matrix scan.
#' Collect treatment and control counts before deciding which pairs pass QC
#' @inheritParams storage_pair_counts
#' @return The count list from `storage_pair_counts()`, in requested pair order,
#'   with optional NT summaries. See that function for which cell numbers to
#'   supply and how stored sparse entries are counted. No QC threshold is applied
#'   here; this step only obtains the measurements.
#' @noRd
collect_pair_qc_metrics <- function(
    store,
    pair_index,
    groups,
    source_cells,
    controls,
    nt_groups = NULL
) {
    storage_pair_counts(
        store,
        pair_index,
        groups,
        source_cells,
        controls,
        nt_groups
    )
}

#' Decide which pairs have enough treatment and control observations
#' @param metrics Count list with integer `n_nonzero_trt` and `n_nonzero_cntrl`
#'   vectors, both in the same pair order.
#' @param qc_options Validated thresholds `n_nonzero_trt_thresh` and
#'   `n_nonzero_cntrl_thresh`, each a single nonmissing integer.
#' @return One logical value per pair: both counts must reach their threshold.
#'   Missing counts follow R's `&` rule, so FALSE & NA is FALSE. Unequal vector
#'   lengths error; empty input gives `logical()`. Calibration sampling uses
#'   this same C++ threshold rule to keep its eligibility decision consistent.
#' @noRd
evaluate_pair_eligibility <- function(metrics, qc_options) {
    evaluate_pair_eligibility_cpp(
        metrics$n_nonzero_trt,
        metrics$n_nonzero_cntrl,
        qc_options$n_nonzero_trt_thresh,
        qc_options$n_nonzero_cntrl_thresh
    )
}
