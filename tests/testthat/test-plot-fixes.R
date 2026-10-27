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

test_that("high gRNA count thresholds do not request negative expansion bins", {
    object <- plot_fix_object()
    before <- serialize(object, NULL)
    set.seed(184L)
    seed <- .Random.seed
    for (threshold in c(16, 40)) {
        panels <- plot_grna_count_distributions(
            object,
            grnas_to_plot = "g",
            threshold = threshold,
            return_indiv_plots = TRUE
        )
        expect_length(panels, 1L)
        expect_identical(
            as.character(panels[[1L]]$data$grna_expressions_bin),
            c("0", "1", "11", "15")
        )
        expect_identical(panels[[1L]]$data$bin_counts, rep(1L, 4L))
    }
    # the exponential bins of the default are unchanged
    default <- plot_grna_count_distributions(
        object,
        grnas_to_plot = "g",
        return_indiv_plots = TRUE
    )
    expect_identical(
        as.character(default[[1L]]$data$grna_expressions_bin),
        c("0", "1", "10", "11")
    )
    expect_identical(default[[1L]]$data$bin_counts, rep(1L, 4L))
    expect_identical(.Random.seed, seed)
    expect_identical(serialize(object, NULL), before)
})

test_that("volcano colors follow the significance calls, including ties", {
    result <- data.frame(
        p_value = c(.01, .02, .02, .03),
        log_2_fold_change = c(-2, 0, .5, 2),
        significant = c(TRUE, TRUE, FALSE, FALSE)
    )
    before <- serialize(result, NULL)
    plot <- make_volcano_plot(result, p_thresh = .02)
    expect_identical(plot$data$reject, result$significant)
    expect_identical(plot$data$p_value, result$p_value)
    expect_identical(plot$data$log_2_fold_change, c(-1.5, 0, .5, 1.5))
    expect_identical(plot$layers[[2L]]$data$yintercept, .02)
    expect_identical(serialize(result, NULL), before)
    colors <- function(plot) ggplot2::ggplot_build(plot)$data[[1L]]$colour
    expect_identical(
        colors(plot),
        c("blueviolet", "blueviolet", "dodgerblue", "dodgerblue")
    )
    # the colors do not depend on which calls are present
    result$significant <- TRUE
    expect_identical(
        colors(make_volcano_plot(result, p_thresh = .03)),
        rep("blueviolet", 4L)
    )
    result$significant <- FALSE
    expect_identical(
        colors(make_volcano_plot(result, p_thresh = NA)),
        rep("dodgerblue", 4L)
    )
})
