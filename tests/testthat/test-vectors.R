# gRNAs on one vector are collapsed into a single unit at import. The gRNAs of
# a vector target different coordinates; the collapsed gRNA target data frame
# must still have one row per vector, or every vector is analyzed once per
# gRNA (Bonferroni p-values multiplied, singleton rows repeated).

vector_target_data_frame <- function() {
    data(grna_target_data_frame_highmoi, package = "sceptre", envir = environment())
    targets <- as.data.frame(grna_target_data_frame_highmoi)
    # two gRNAs per vector, each vector within one target
    targets$vector_id <- paste0(
        "vector_",
        targets$grna_target,
        "_",
        stats::ave(
            seq_len(nrow(targets)),
            targets$grna_target,
            FUN = function(i) (seq_along(i) + 1L) %/% 2L
        )
    )
    targets
}


test_that("vectors with target coordinates collapse to one row each", {
    targets <- vector_target_data_frame()
    collapsed <- collapse_grna_target_data_frame(targets)
    expect_identical(sort(collapsed$grna_id), sort(unique(targets$vector_id)))
    expect_false(anyDuplicated(collapsed$grna_id) > 0L)
    # the region of each vector spans those of its gRNAs
    expected_start <- tapply(targets$start, targets$vector_id, min)
    expected_end <- tapply(targets$end, targets$vector_id, max)
    expect_equal(
        collapsed$start,
        as.vector(expected_start[collapsed$grna_id])
    )
    expect_equal(collapsed$end, as.vector(expected_end[collapsed$grna_id]))
    expect_s3_class(collapsed, "data.frame")
    expect_false(inherits(collapsed, "tbl_df"))

    # a vector whose gRNAs target different chromosomes is an error
    on_two_chromosomes <- targets
    on_two_chromosomes$chr[2] <- "chrX_other"
    expect_error(
        collapse_grna_target_data_frame(on_two_chromosomes),
        "one chromosome"
    )

    # columns that describe the gRNAs rather than the vector are dropped
    annotated <- targets
    annotated$spacer <- paste0("spacer_", seq_len(nrow(annotated)))
    annotated$library <- "library_1"
    expect_message(
        collapsed_annotated <- collapse_grna_target_data_frame(annotated),
        "`spacer`"
    )
    expect_false("spacer" %in% colnames(collapsed_annotated))
    expect_true("library" %in% colnames(collapsed_annotated))
    expect_identical(nrow(collapsed_annotated), nrow(collapsed))

    # a factor vector_id gives the collapsed gRNA matrix its row names
    data(highmoi_example_data, package = "sceptre", envir = environment())
    with_factor <- targets[, c("grna_id", "grna_target", "vector_id")]
    with_factor$vector_id <- factor(with_factor$vector_id)
    collapsed_matrix <- collapse_grna_matrix(
        highmoi_example_data$grna_matrix,
        with_factor
    )
    expect_identical(
        rownames(collapsed_matrix),
        as.character(unique(with_factor$vector_id))
    )
})


test_that("target coordinates do not change the analysis of vectors", {
    data(highmoi_example_data, package = "sceptre", envir = environment())
    targets <- vector_target_data_frame()
    without_coordinates <- targets[, c("grna_id", "grna_target", "vector_id")]
    analyze <- function(grna_target_data_frame) {
        sceptre_object <- import_data(
            response_matrix = highmoi_example_data$response_matrix,
            grna_matrix = highmoi_example_data$grna_matrix,
            grna_target_data_frame = grna_target_data_frame,
            moi = "high",
            extra_covariates = highmoi_example_data$extra_covariates
        )
        pairs <- data.frame(
            grna_target = rep(c("candidate_enh_1", "candidate_enh_2"), each = 3),
            response_id = rep(
                rownames(highmoi_example_data$response_matrix)[1:3],
                2
            )
        )
        suppressMessages(
            sceptre_object |>
                set_analysis_parameters(
                    discovery_pairs = pairs,
                    side = "left",
                    grna_integration_strategy = "bonferroni",
                    resampling_mechanism = "permutations",
                    seed = 4
                ) |>
                assign_grnas(method = "thresholding", print_progress = FALSE) |>
                run_qc(n_nonzero_trt_thresh = 0L, n_nonzero_cntrl_thresh = 0L) |>
                run_discovery_analysis(print_progress = FALSE)
        )
    }
    with_coordinates <- analyze(targets)
    without <- analyze(without_coordinates)
    expect_false(anyDuplicated(with_coordinates@grna_target_data_frame$grna_id) > 0L)
    expect_equal(
        get_result(with_coordinates, "run_discovery_analysis"),
        get_result(without, "run_discovery_analysis")
    )
    # the gRNA groups, and so every analysis (including the default
    # calibration group size, which counts the rows of each target), are
    # those without coordinates
    columns <- c("grna_id", "grna_target")
    expect_equal(
        as.data.frame(with_coordinates@grna_target_data_frame)[, columns],
        as.data.frame(without@grna_target_data_frame)[, columns]
    )
    expect_identical(
        compute_calibration_group_size(with_coordinates@grna_target_data_frame),
        compute_calibration_group_size(without@grna_target_data_frame)
    )
})
