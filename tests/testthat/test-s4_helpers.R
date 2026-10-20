test_that("display and generic plot derive state from completed components", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$low_complement_perm
    imported <- refactor_capture(refactor_import(inputs, case))$value
    expect_identical(slotNames(imported),
        c("data", "spec", "preparation", "fits", "results", "schema_version"))
    expect_output(show(imported), "cells")
    expect_output(print(imported), "Analysis status")
    expect_error(plot(imported), "no generic plot", ignore.case = TRUE)
    object <- refactor_capture(refactor_prepare(inputs, case))$value
    before <- serialize(object, NULL)
    testthat::local_mocked_bindings(
        plot_run_qc = function(...) "qc",
        plot_run_calibration_check = function(...) "calibration",
        plot_run_discovery_analysis = function(...) "discovery", .package = "sceptre")
    expect_identical(plot(object), "qc")
    calibrated <- refactor_capture(run_calibration_check(object,
        n_calibration_pairs = 2L, calibration_group_size = 1L,
        print_progress = FALSE))$value
    expect_identical(plot(calibrated), "calibration")
    discovered <- refactor_capture(run_discovery_analysis(calibrated,
        print_progress = FALSE))$value
    expect_identical(plot(discovered), "discovery")
    expect_output(print(discovered), "Summary of results")
    changed <- suppressMessages(assign_grnas(discovered, method = "thresholding",
        threshold = 45, print_progress = FALSE))
    testthat::local_mocked_bindings(plot_assign_grnas = function(...) "assignment", .package = "sceptre")
    expect_identical(plot(changed), "assignment")
    expect_identical(serialize(object, NULL), before)
})
