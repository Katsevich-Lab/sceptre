# The Screen/Analysis objects and their files. These tests cover sceptre's own
# mapping; SingleCellExperiment, MuData, anndataR and arrow test their formats.

analysis_workflow <- function() {
    data(lowmoi_example_data, package = "sceptre", envir = environment())
    object <- import_data(
        response_matrix = lowmoi_example_data$response_matrix,
        grna_matrix = lowmoi_example_data$grna_matrix,
        extra_covariates = lowmoi_example_data$extra_covariates,
        grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
        moi = "low"
    )
    positive <- construct_positive_control_pairs(object)
    discovery <- construct_trans_pairs(
        object,
        positive_control_pairs = positive,
        pairs_to_exclude = "pc_pairs"
    )[1:30, ]
    steps <- list(
        parameters = function(o) {
            set_analysis_parameters(
                o,
                discovery_pairs = discovery,
                positive_control_pairs = positive
            )
        },
        assign = function(o) assign_grnas(o, method = "thresholding"),
        qc = function(o) run_qc(o),
        discovery = function(o) {
            run_discovery_analysis(o, print_progress = FALSE)
        }
    )
    list(start = object, steps = steps)
}

run_steps <- function(object, steps) {
    for (step in steps) {
        object <- step(object)
    }
    object
}

test_that("every sceptre_object slot belongs to exactly one part", {
    mapped <- unname(unlist(analysis_layout))
    expect_setequal(mapped, methods::slotNames("sceptre_object"))
    expect_false(anyDuplicated(mapped) > 0L)
})

test_that("converting to an Analysis and back loses nothing at any stage", {
    workflow <- analysis_workflow()
    object <- workflow$start
    expect_identical(as_sceptre_object(as_analysis(object)), object)
    for (step in workflow$steps) {
        object <- suppressMessages(step(object))
        expect_identical(as_sceptre_object(as_analysis(object)), object)
    }
    screen <- as_analysis(object)@screen
    expect_s4_class(screen, "SingleCellExperiment")
    expect_identical(SingleCellExperiment::altExpNames(screen), "gRNA")
})

test_that("an analysis saved to files and read back finishes with the same results", {
    skip_if_not_installed("anndataR")
    skip_if_not_installed("arrow")
    skip_if_not_installed("MuData")
    workflow <- analysis_workflow()
    steps <- workflow$steps
    direct <- suppressMessages(run_steps(workflow$start, steps))
    # save after QC, then run the discovery analysis on what was read back
    before <- suppressMessages(run_steps(workflow$start, steps[1:3]))
    directory <- withr::local_tempdir()
    write_analysis(as_analysis(before), directory)
    expect_true(all(file.exists(file.path(
        directory,
        c("screen.h5mu", "state.h5ad")
    ))))
    read_back <- as_sceptre_object(read_analysis(directory))
    resumed <- suppressMessages(steps$discovery(read_back))
    expect_identical(resumed@discovery_result, direct@discovery_result)
    # results and fits are read back from Parquet with their values intact
    write_analysis(as_analysis(direct), directory)
    final <- read_analysis(directory)
    expect_equal(
        as.data.frame(final@results$discovery),
        as.data.frame(direct@discovery_result),
        ignore_attr = TRUE
    )
    expect_equal(final@fits, direct@response_precomputations)
})

test_that("missing gRNA locations survive the screen file", {
    skip_if_not_installed("MuData")
    data(highmoi_example_data, package = "sceptre", envir = environment())
    data(
        grna_target_data_frame_highmoi,
        package = "sceptre",
        envir = environment()
    )
    object <- import_data(
        response_matrix = highmoi_example_data$response_matrix,
        grna_matrix = highmoi_example_data$grna_matrix,
        extra_covariates = highmoi_example_data$extra_covariates,
        grna_target_data_frame = grna_target_data_frame_highmoi,
        moi = "high"
    )
    expect_true(anyNA(object@grna_target_data_frame$chr))
    directory <- withr::local_tempdir()
    write_analysis(as_analysis(object), directory)
    read_back <- as_sceptre_object(read_analysis(directory))
    expect_identical(
        read_back@grna_target_data_frame,
        object@grna_target_data_frame
    )
    expect_identical(
        read_back@covariate_data_frame,
        object@covariate_data_frame
    )
})
