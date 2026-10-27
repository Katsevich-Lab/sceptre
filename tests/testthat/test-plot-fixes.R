plot_fix_object <- function() {
    object <- methods::new("sceptre_object")
    object <- set_grna_matrix(
        object,
        matrix(c(0, 1, 11, 15), 1L, dimnames = list("g", NULL))
    )
    object <- set_response_matrix(object, matrix(1, 1L, 4L))
    object@low_moi <- FALSE
    object@cell_removal_metrics <- rep(0L, 6L)
    object@functs_called <- c(run_qc = TRUE)
    object@n_nonzero_cntrl_thresh <- object@n_nonzero_trt_thresh <- 7L
    object@discovery_pairs_with_info <- data.frame(
        response_id = paste0("r", 1:6),
        grna_group = "g",
        pass_qc = rep(c(TRUE, FALSE), 3L),
        n_nonzero_trt = 1:6,
        n_nonzero_cntrl = 6:1
    )
    object
}

test_that("plot_run_qc() forwards its pairwise plot options", {
    object <- plot_fix_object()
    before <- serialize(object, NULL)
    set.seed(829L)
    seed <- .Random.seed
    panels <- plot_run_qc(
        object,
        downsample_pairs = 2L,
        point_size = 3,
        transparency = .2,
        return_indiv_plots = TRUE
    )
    expect_length(panels, 2L)
    expected_rows <- withr::with_seed(4, sample.int(6L, 2L))
    expect_identical(
        panels[[2L]]$data$response_id,
        object@discovery_pairs_with_info$response_id[expected_rows]
    )
    expect_identical(panels[[2L]]$layers[[1L]]$aes_params$size, 3)
    expect_identical(panels[[2L]]$layers[[1L]]$aes_params$alpha, .2)
    # defaults are unchanged
    default <- plot_run_qc(object, return_indiv_plots = TRUE)[[2L]]
    direct <- plot_pairwise_qc(object)
    expect_identical(default$data, direct$data)
    expect_identical(
        default$layers[[1L]]$aes_params,
        direct$layers[[1L]]$aes_params
    )
    expect_identical(nrow(default$data), 6L)
    object@discovery_pairs_with_info <- object@discovery_pairs_with_info[
        FALSE,
    ]
    single <- plot_run_qc(object, return_indiv_plots = TRUE)
    expect_s3_class(single, "ggplot")
    expect_identical(.Random.seed, seed)
    object@discovery_pairs_with_info <- unserialize(
        before
    )@discovery_pairs_with_info
    expect_identical(serialize(object, NULL), before)
})
