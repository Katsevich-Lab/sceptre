# Measurement and policy are independent of the S4 facade. The storage adapter
# retains its fused NT summary scan; policy consumes only aligned count vectors.
collect_pair_qc_metrics <- function(store, pair_index, groups, source_cells,
                                    controls, nt_groups = NULL, nuclear = FALSE) {
    storage_pair_counts(store, pair_index, groups, source_cells,
                        controls, nt_groups, nuclear)
}

evaluate_pair_eligibility <- function(metrics, qc_options) {
    evaluate_pair_eligibility_cpp(
        metrics$n_nonzero_trt, metrics$n_nonzero_cntrl,
        qc_options$n_nonzero_trt_thresh, qc_options$n_nonzero_cntrl_thresh
    )
}
