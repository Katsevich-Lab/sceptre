# A sceptre_object written in one R session must be usable in another in which
# only sceptre is loaded. Its data frame slots may hold tibbles (the gRNA target
# data frame of the example data) and data.tables (pairs and results), which
# S4 recognizes as data frames only once tibble and data.table are loaded.
# These tests read objects in a separate R process; the tests run in a session
# in which those packages are already loaded could not detect the problem.

analyzed_example_object <- function() {
    data(lowmoi_example_data, package = "sceptre", envir = environment())
    sceptre_object <- import_data(
        response_matrix = lowmoi_example_data$response_matrix,
        grna_matrix = lowmoi_example_data$grna_matrix,
        grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
        extra_covariates = lowmoi_example_data$extra_covariates,
        moi = "low"
    )
    suppressMessages(
        sceptre_object |>
            set_analysis_parameters(
                discovery_pairs = construct_trans_pairs(sceptre_object),
                seed = 4
            ) |>
            assign_grnas(print_progress = FALSE) |>
            run_qc() |>
            run_calibration_check(
                n_calibration_pairs = 20,
                print_progress = FALSE
            )
    )
}


test_that("read_sceptre_object() works in a fresh R session", {
    skip_if_not_installed("callr")
    sceptre_object <- analyzed_example_object()
    expect_s3_class(sceptre_object@grna_target_data_frame, "tbl_df")
    expect_s3_class(sceptre_object@calibration_result, "data.table")
    directory <- withr::local_tempdir()
    write_sceptre_object(sceptre_object, directory)

    out <- run_in_fresh_session(
        function(directory) {
            result <- function(backed) {
                x <- sceptre::read_sceptre_object(directory, backed = backed)
                methods::validObject(x)
                sceptre::get_result(x, "run_calibration_check")
            }
            list(memory = result(FALSE), backed = result(TRUE))
        },
        directory = directory
    )
    # the test is meaningful only if the other process had not loaded them
    skip_if(
        any(c("data.table", "tibble") %in% out$loaded_before),
        "library(sceptre) now loads data.table or tibble"
    )
    expected <- get_result(sceptre_object, "run_calibration_check")
    expect_equal(out$value$memory, expected)
    expect_equal(out$value$backed, expected)
})


test_that("updateObject() works in a fresh R session", {
    skip_if_not_installed("callr")
    sceptre_object <- analyzed_example_object()
    file_old <- withr::local_tempfile(fileext = ".rds")
    file_current <- withr::local_tempfile(fileext = ".rds")
    saveRDS(old_layout_object(sceptre_object), file_old)
    saveRDS(sceptre_object, file_current)

    out <- run_in_fresh_session(
        function(file_old, file_current) {
            update <- function(file) {
                x <- BiocGenerics::updateObject(readRDS(file))
                sceptre::get_result(x, "run_calibration_check")
            }
            list(old = update(file_old), current = update(file_current))
        },
        file_old = file_old,
        file_current = file_current
    )
    skip_if(
        any(c("data.table", "tibble") %in% out$loaded_before),
        "library(sceptre) now loads data.table or tibble"
    )
    expected <- get_result(sceptre_object, "run_calibration_check")
    expect_equal(out$value$old, expected)
    expect_equal(out$value$current, expected)
})


test_that("updateObject() converts an object saved by sceptre 0.99.1", {
    # fixtures/sceptre_0.99.1_lowmoi.rds was saved by sceptre 0.99.1 (its gRNA
    # target data frame is a tibble, its calibration result a data.table):
    # the low-MOI example data, analyzed through run_calibration_check()
    # with n_calibration_pairs = 20
    skip_if_not_installed("callr")
    file <- test_path("fixtures", "sceptre_0.99.1_lowmoi.rds")
    out <- run_in_fresh_session(
        function(file) {
            x <- BiocGenerics::updateObject(readRDS(file))
            methods::validObject(x)
            list(
                version = as.character(x@version),
                seed = x@seed,
                calibration = sceptre::get_result(x, "run_calibration_check"),
                n_cells = ncol(sceptre::get_response_matrix(x))
            )
        },
        file = normalizePath(file)
    )
    old <- readRDS(file)
    expect_identical(out$value$version, "0.99.1")
    expect_identical(out$value$seed, 4L)
    expect_identical(out$value$n_cells, ncol(attr(old, "response_matrix")[[1]]))
    expect_equal(
        as.data.frame(out$value$calibration),
        as.data.frame(attr(old, "calibration_result"))
    )
})
