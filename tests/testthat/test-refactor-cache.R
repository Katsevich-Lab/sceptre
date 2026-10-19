test_that("discovery calibration discovery preserves valid fit views", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    cases <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))
    for (case in cases[c("low_nt_perm", "low_complement_perm", "high_complement_crt")]) {
        object <- refactor_prepare(inputs, case)
        discovery <- suppressMessages(sceptre::run_discovery_analysis(
            object, print_progress = FALSE, output_amount = 2L))
        calibrated <- suppressMessages(sceptre::run_calibration_check(
            discovery, n_calibration_pairs = 4L, calibration_group_size = 2L,
            print_progress = FALSE, output_amount = 2L))
        cold <- calibrated
        cold@fits@response_precomputations <- list()
        warm <- suppressMessages(sceptre::run_discovery_analysis(
            calibrated, print_progress = FALSE, output_amount = 2L))
        cold <- suppressMessages(sceptre::run_discovery_analysis(
            cold, print_progress = FALSE, output_amount = 2L))
        expect_equal(warm@results@discovery_result, cold@results@discovery_result, tolerance = 1e-8)
        expect_identical(object@fits@response_precomputations, list())
    }
})

test_that("changing the control view discards incompatible calibration fits", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$low_nt_perm
    object <- refactor_prepare(inputs, case)
    object <- suppressMessages(sceptre::run_calibration_check(object,
        n_calibration_pairs = 4L, calibration_group_size = 2L,
        print_progress = FALSE, output_amount = 2L))
    expect_gt(length(object@fits@response_precomputations), 0L)
    changed <- sceptre::set_analysis_parameters(object,
        discovery_pairs = inputs$discovery_pairs,
        positive_control_pairs = inputs$positive_control_pairs,
        formula_object = ~ batch + z, control_group = "complement",
        resampling_mechanism = "permutations")
    expect_length(changed@fits@response_precomputations, 0L)
    changed <- sceptre::assign_grnas(changed, method = "thresholding",
        threshold = 40, print_progress = FALSE)
    changed <- sceptre::run_qc(changed, response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1), n_nonzero_trt_thresh = 7L,
        n_nonzero_cntrl_thresh = 7L, additional_cells_to_remove = inputs$remove_cells)
    cold <- changed
    cold@fits@response_precomputations <- list()
    warm <- suppressMessages(sceptre::run_discovery_analysis(changed,
        print_progress = FALSE, output_amount = 2L))
    cold <- suppressMessages(sceptre::run_discovery_analysis(cold,
        print_progress = FALSE, output_amount = 2L))
    expect_equal(warm@results@discovery_result, cold@results@discovery_result, tolerance = 1e-8)
})

test_that("actual design contrasts rather than formula text govern cache reuse", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$low_complement_perm
    object <- refactor_prepare(inputs, case)
    object <- suppressMessages(sceptre::run_discovery_analysis(object,
        print_progress = FALSE, output_amount = 2L))
    expect_gt(length(object@fits@response_precomputations), 0L)
    before <- serialize(object, NULL)
    changed <- withr::with_options(list(contrasts = c("contr.sum", "contr.poly")),
        sceptre::set_analysis_parameters(object,
            discovery_pairs = inputs$discovery_pairs,
            formula_object = ~ batch + z, control_group = "complement",
            resampling_mechanism = "permutations"))
    expect_false(identical(changed@preparation@covariate_matrix, object@preparation@covariate_matrix))
    expect_length(changed@fits@response_precomputations, 0L)
    expect_identical(serialize(object, NULL), before)
})
