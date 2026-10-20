test_that("workflow status follows completed components and preserves sibling results", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$low_complement_perm
    object <- refactor_capture(refactor_prepare(inputs, case))$value
    original <- serialize(object, NULL)
    status <- workflow_status(object)
    expect_identical(names(status), c("import_data", "set_analysis_parameters", "assign_grnas", "run_qc",
        "run_calibration_check", "run_power_check", "run_discovery_analysis"))
    expect_identical(unname(status), c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))
    discovery <- refactor_capture(run_discovery_analysis(object, print_progress = FALSE))$value
    expect_true(workflow_status(discovery)[["run_discovery_analysis"]])
    calibrated <- refactor_capture(run_calibration_check(discovery,
        n_calibration_pairs = 2L, calibration_group_size = 1L,
        print_progress = FALSE))$value
    expect_true(all(workflow_status(calibrated)[c("run_discovery_analysis", "run_calibration_check")]))
    expect_false(workflow_status(calibrated)[["run_power_check"]])
    changed <- suppressMessages(assign_grnas(calibrated, method = "thresholding",
        threshold = 45, print_progress = FALSE))
    expect_true(workflow_status(changed)[["assign_grnas"]])
    expect_false(any(workflow_status(changed)[c("run_qc", "run_discovery_analysis", "run_calibration_check")]))
    expect_identical(serialize(object, NULL), original)
})
