result_test_object <- function() {
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    case <- readRDS(test_path("fixtures", "refactor", "cases.rds"))$low_complement_perm
    inputs$discovery_pairs <- inputs$discovery_pairs[3L, , drop = FALSE]
    object <- suppressMessages(refactor_prepare(inputs, case))
    object@spec@multiple_testing_alpha <- 0.1
    object
}

result_test_commit <- function(object, p = 0.2) {
    commit_result_record(object, "run_discovery_analysis",
        data.frame(response_id = "r1", grna_target = "A", p_value = p,
            significant = p < object@spec@multiple_testing_alpha),
        object@preparation@discovery_pairs_with_info, 1,
        plan = make_run_plan(object, "run_discovery_analysis",
            budgets = c(B1 = 2L, B2 = 3L, B3 = 5L)))
}

test_that("completed results freeze effective budgets and preserve historical tables", {
    object <- result_test_commit(result_test_object())
    expect_true(result_record_is_current(object, "run_discovery_analysis"))
    expect_identical(get_result(object, "run_discovery_analysis")$p_value, 0.2)
    alternate <- make_run_plan(object, "run_discovery_analysis",
        budgets = c(B1 = 2L, B2 = 3L, B3 = 99L))
    expect_identical(alternate@B3, 99L)
    expect_true(result_record_is_current(object, "run_discovery_analysis"))
    expect_identical(object@results@records$run_discovery_analysis$provenance$plan@B3, 5L)
    object <- invalidate_result_records(object, "design changed")
    expect_error(get_result(object, "run_discovery_analysis"), "No completed current result")
    expect_warning(old <- get_result(object, "run_discovery_analysis", historical = TRUE), "historical")
    expect_identical(old$p_value, 0.2)
    expect_identical(attr(old, "sceptre_result_record")$invalidated_by, "design changed")
    object <- result_test_commit(object, 0.4)
    expect_identical(object@results@records$run_discovery_analysis$history[[1]]$table$p_value, 0.2)
    expect_identical(get_result(object, "run_discovery_analysis")$p_value, 0.4)
})

test_that("result guards detect data, ordered cell, design and family changes", {
    original <- result_test_commit(result_test_object())
    for (mutation in c("cells", "design", "data", "pairs", "alpha")) {
        object <- original
        if (mutation == "cells") object@preparation@analysis_cells <- rev(object@preparation@analysis_cells)
        if (mutation == "design") object@preparation@design_matrix[3, 2] <- 4
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
    prepared <- invalidate_result_records(object, "analysis prepared for rerun")
    expect_false(workflow_status(prepared)[["run_discovery_analysis"]])
    expect_error(get_result(prepared, "run_discovery_analysis"), "No completed current result")
    expect_true(workflow_status(object)[["run_discovery_analysis"]])
    shard <- object
    shard@results@records$run_discovery_analysis$scope <- "raw_shard"
    expect_error(get_result(shard, "run_discovery_analysis"), "No completed current result")
})

test_that("zero discovery eligibility does not prevent qualifying NT calibration proposals", {
    # Nine NT guides choose four = 126 groups: exercise the proposal sampler,
    # rather than the small-combination enumeration branch.
    guide_ids <- c(paste0("nt", 1:9), "target-guide")
    guide <- matrix(0, 10L, 200L,
        dimnames = list(guide_ids, paste0("c", 1:200)))
    guide[cbind(rep(1:10, each = 20L), 1:200)] <- 50
    object <- suppressMessages(import_data(
        response_matrix = matrix(1, 1L, 200L,
            dimnames = list("r1", colnames(guide))), grna_matrix = guide,
        grna_target_data_frame = data.frame(grna_id = guide_ids,
            grna_target = c(rep("non-targeting", 9L), "T")), moi = "low"))
    object <- set_analysis_parameters(object,
        discovery_pairs = data.frame(response_id = "r1", grna_target = "T"),
        formula_object = ~ 1, control_group = "complement",
        resampling_approximation = "no_approximation")
    object <- suppressMessages(assign_grnas(object, method = "thresholding",
        threshold = 40, print_progress = FALSE))
    object <- suppressMessages(run_qc(object, response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1), n_nonzero_trt_thresh = 50L,
        n_nonzero_cntrl_thresh = 50L))
    expect_false(any(object@preparation@discovery_pairs_with_info$pass_qc))
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
    expect_false(workflow_status(prepared$object)[["run_calibration_check"]])
    expect_null(prepared$object@results@records$run_calibration_check)
    expect_gt(nrow(prepared$plan@pairs), 0L)
    before <- serialize(prepared, NULL)
    testthat::local_mocked_bindings(analyze_association_pairs = function(...) stop("worker failure"))
    expect_error(run_calibration_check_pt_2(prepared$object, prepared$plan, print_progress = FALSE), "worker failure")
    expect_identical(serialize(prepared, NULL), before)
    expect_false(workflow_status(prepared$object)[["run_calibration_check"]])
    expect_null(prepared$object@results@records$run_calibration_check)
})

test_that("calibration working budget does not leak into power or discovery", {
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    case <- readRDS(test_path("fixtures", "refactor", "cases.rds"))$low_complement_perm
    object <- suppressMessages(refactor_prepare(inputs, case))
    expected <- make_run_plan(object, "run_discovery_analysis")@B3
    prepared <- suppressMessages(run_calibration_check_pt_1(object,
        n_calibration_pairs = 1L, calibration_group_size = 1L))
    expect_false(identical(prepared$plan@B3, expected))
    object <- prepared$object
    observed <- list()
    testthat::local_mocked_bindings(analyze_association_pairs = function(sceptre_object, ..., plan) {
        observed[[length(observed) + 1L]] <<- plan@B3
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
    expect_false(workflow_status(object)[["run_discovery_analysis"]])
})


test_that("result records preserve the public result table's data.table class", {
    object <- result_test_object()
    result <- data.table::data.table(response_id = "r1", grna_target = "A", p_value = 0.2,
        significant = FALSE)
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
    frozen <- object@results@records$run_discovery_analysis$provenance
    expect_identical(frozen$spec$grna_assignment_method, "thresholding")
    expect_identical(frozen$spec$grna_assignment_hyperparameters$threshold, 40)
    expect_identical(frozen$spec$n_nonzero_trt_thresh, 7L)
    expect_identical(frozen$response_fit_method, "poisson-glm.fit-then-dispersion-v1")
})
