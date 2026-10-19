result_test_object <- function() {
    object <- methods::new("sceptre_object")
    object@functs_called <- workflow_initial_state()
    object@functs_called[c("set_analysis_parameters", "assign_grnas", "run_qc")] <- TRUE
    object@data@identity <- "test-content-v1"
    object@preparation@cells_in_use <- c(2L, 7L, 9L)
    object@preparation@covariate_matrix <- matrix(c(1, 1, 1, 0, 1, 2), 3L,
        dimnames = list(c("c2", "c7", "c9"), c("intercept", "z")))
    object@preparation@discovery_pairs_with_info <- data.frame(response_id = "r1", grna_group = "g1", pass_qc = TRUE)
    object@spec@B1 <- 2L
    object@spec@B2 <- 3L
    object@spec@B3 <- 5L
    object@spec@multiple_testing_alpha <- 0.1
    object@spec@multiple_testing_method <- "BH"
    object@nf_pipeline <- FALSE
    object
}

result_test_commit <- function(object, p = 0.2) {
    commit_result_record(object, "run_discovery_analysis",
        data.frame(response_id = "r1", grna_target = "g1", p_value = p),
        object@preparation@discovery_pairs_with_info, 1)
}

test_that("completed results freeze effective budgets and preserve historical tables", {
    object <- result_test_commit(result_test_object())
    expect_true(result_record_is_current(object, "run_discovery_analysis"))
    expect_identical(get_result(object, "run_discovery_analysis")$p_value, 0.2)
    object@spec@B3 <- 99L
    expect_true(result_record_is_current(object, "run_discovery_analysis"))
    expect_identical(object@results@records$run_discovery_analysis$provenance$effective$B3, 5L)
    object <- invalidate_result_records(object, "design changed")
    expect_error(get_result(object, "run_discovery_analysis"), "No completed current result")
    expect_warning(old <- get_result(object, "run_discovery_analysis", historical = TRUE), "historical")
    expect_identical(old$p_value, 0.2)
    expect_identical(attr(old, "sceptre_result_record")$invalidation_reason, "design changed")
    object <- result_test_commit(object, 0.4)
    expect_identical(object@results@records$run_discovery_analysis$history[[1]]$result$p_value, 0.2)
    expect_identical(get_result(object, "run_discovery_analysis")$p_value, 0.4)
})

test_that("result guards detect data, ordered cell, design and family changes", {
    original <- result_test_commit(result_test_object())
    for (mutation in c("cells", "design", "data", "pairs", "alpha")) {
        object <- original
        if (mutation == "cells") object@preparation@cells_in_use <- rev(object@preparation@cells_in_use)
        if (mutation == "design") object@preparation@covariate_matrix[3, 2] <- 4
        if (mutation == "data") object@data@identity <- "changed-content"
        if (mutation == "pairs") object@preparation@discovery_pairs_with_info$response_id <- "r2"
        if (mutation == "alpha") object@spec@multiple_testing_alpha <- 0.05
        expect_false(result_record_is_current(object, "run_discovery_analysis"), info = mutation)
        expect_error(get_result(object, "run_discovery_analysis"), "No completed current result")
        expect_equal(current_result(object, "run_discovery_analysis", required = FALSE), data.frame())
    }
})

test_that("prepared analysis and raw shards cannot masquerade as completed reporting results", {
    object <- result_test_commit(result_test_object())
    prepared <- begin_analysis_stage(object, "run_discovery_analysis")
    expect_false(prepared@functs_called[["run_discovery_analysis"]])
    expect_error(get_result(prepared, "run_discovery_analysis"), "No completed current result")
    expect_true(object@functs_called[["run_discovery_analysis"]])
    object@nf_pipeline <- TRUE
    shard <- result_test_commit(object)
    expect_error(get_result(shard, "run_discovery_analysis"), "raw shards")
})

test_that("zero discovery eligibility does not prevent qualifying NT calibration proposals", {
    # Nine NT guides choose four = 126 groups: exercise the proposal sampler,
    # rather than the small-combination enumeration branch.
    object <- result_test_object()
    object@data@low_moi <- TRUE
    object@data@response_matrix <- list(Matrix::Matrix(matrix(1, 1, 200,
        dimnames = list("r1", paste0("c", 1:200))), sparse = TRUE))
    object@preparation@grna_assignments <- list(indiv_nt_grna_idxs =
        stats::setNames(lapply(0:8, function(i) as.integer(i * 20 + 1:20)), paste0("nt", 1:9)))
    object@preparation@cells_in_use <- 1:200
    object@preparation@M_matrix <- matrix(20L, 9, 1)
    object@preparation@n_nonzero_tot_vector <- 200L
    object@preparation@n_discovery_pairs <- 1L
    object@preparation@discovery_pairs_with_info <- data.frame(pass_qc = FALSE)
    object@spec@n_nonzero_trt_thresh <- 50L
    object@spec@n_nonzero_cntrl_thresh <- 50L
    pairs <- construct_negative_control_pairs(object, 4L, 4L)
    expect_equal(nrow(pairs), 4L)
    expect_true(all(pairs$pass_qc))
    expect_true(all(pairs$n_nonzero_trt == 80L))
    expect_true(all(pairs$n_nonzero_cntrl == 120L))
})

test_that("public calibration preparation stays incomplete and failed execution commits nothing", {
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    case <- readRDS(test_path("fixtures", "refactor", "cases.rds"))$low_complement_perm
    object <- suppressMessages(refactor_prepare(inputs, case))
    prepared <- suppressMessages(run_calibration_check_pt_1(object,
        n_calibration_pairs = 2L, calibration_group_size = 1L))
    expect_false(prepared@functs_called[["run_calibration_check"]])
    expect_null(prepared@results@records$run_calibration_check)
    expect_gt(nrow(prepared@spec@negative_control_pairs), 0L)
    before <- serialize(prepared, NULL)
    testthat::local_mocked_bindings(analyze_association_pairs = function(...) stop("worker failure"))
    expect_error(run_calibration_check_pt_2(prepared, print_progress = FALSE), "worker failure")
    expect_identical(serialize(prepared, NULL), before)
    expect_false(prepared@functs_called[["run_calibration_check"]])
    expect_null(prepared@results@records$run_calibration_check)
})

test_that("calibration working budget does not leak into power or discovery", {
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    case <- readRDS(test_path("fixtures", "refactor", "cases.rds"))$low_complement_perm
    object <- suppressMessages(refactor_prepare(inputs, case))
    expected <- object@spec@B3
    object <- suppressMessages(run_calibration_check_pt_1(object,
        n_calibration_pairs = 1L, calibration_group_size = 1L))
    expect_false(identical(object@spec@B3, expected))
    observed <- list()
    testthat::local_mocked_bindings(analyze_association_pairs = function(sceptre_object, ...) {
        observed[[length(observed) + 1L]] <<- sceptre_object@spec@B3
        stop("budget observed")
    })
    expect_error(suppressMessages(run_power_check(object, print_progress = FALSE)), "budget observed")
    expect_error(suppressMessages(run_discovery_analysis(object, print_progress = FALSE)), "budget observed")
    expect_identical(observed, list(expected, expected))
})


test_that("invalid result proposals cannot create completed records", {
    object <- result_test_object()
    before <- serialize(object, NULL)
    for (proposal in list(data.frame(), data.frame(response_id = "r1"), list(p_value = 0.1))) {
        expect_error(commit_result_record(object, "run_discovery_analysis", proposal,
            data.frame(), 1), "completed analysis requires")
    }
    expect_identical(serialize(object, NULL), before)
    expect_false(object@functs_called[["run_discovery_analysis"]])
})


test_that("result records preserve the public result table's data.table class", {
    object <- result_test_object()
    result <- data.table::data.table(response_id = "r1", p_value = 0.2)
    object <- commit_result_record(object, "run_discovery_analysis", result,
        object@preparation@discovery_pairs_with_info, 1)
    expect_s3_class(get_result(object, "run_discovery_analysis"), "data.table")
    expect_true(validObject(object@results))
})


test_that("run provenance freezes assignment and QC settings separately from working choices", {
    object <- result_test_object()
    object@spec@grna_assignment_method <- "thresholding"
    object@spec@grna_assignment_hyperparameters <- list(threshold = 40)
    object@spec@n_nonzero_trt_thresh <- 7L
    object <- result_test_commit(object)
    object@spec@grna_assignment_hyperparameters <- list(threshold = 41)
    object@spec@n_nonzero_trt_thresh <- 10L
    effective <- object@results@records$run_discovery_analysis$provenance$effective
    expect_identical(effective$assignment$method, "thresholding")
    expect_identical(effective$assignment$hyperparameters$threshold, 40)
    expect_identical(effective$qc$treatment_nonzero, 7L)
    expect_identical(effective$response_fit_method, "poisson-glm.fit-then-dispersion-v1")
})
