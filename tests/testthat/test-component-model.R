component_workflow <- function(name = "high_complement_crt") {
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    case <- readRDS(test_path("fixtures", "refactor", "cases.rds"))[[name]]
    suppressMessages(refactor_prepare(inputs, case))
}

test_that("component validity rejects inconsistent cross-object coordinates", {
    object <- component_workflow()
    expect_true(validObject(object))
    expect_s4_class(object@data, "ScreenData")
    expect_s4_class(object@spec, "AnalysisSpec")
    expect_s4_class(object@preparation, "PreparationRecord")
    bad <- object
    bad@preparation@cells_in_use <- c(1L, 1L)
    expect_error(validObject(bad), "unique source-cell")
    bad <- object
    bad@data@grna_matrix <- list(get_grna_matrix(bad)[, -1L])
    expect_error(validObject(bad), "aligned cells")
    bad <- object
    bad@preparation@initial_grna_assignment_list[[1]] <- 0L
    expect_error(validObject(bad), "canonical source-cell")
    expect_error(new("CellView", fit_to_analysis = 1:2,
        source_cells = c(7L, 9L), trt_fit_pos = 3L, n_trt = 1L), "inconsistent")
})

test_that("count content replacement invalidates stages while preserving historical tables", {
    object <- component_workflow()
    object <- suppressMessages(run_discovery_analysis(object, print_progress = FALSE))
    historical <- get_result(object, "run_discovery_analysis")
    counts <- get_response_matrix(object)
    counts@x[1] <- counts@x[1] + 1
    changed <- set_response_matrix(object, counts)
    expect_false(identical(changed@data@identity, object@data@identity))
    expect_length(changed@fits@response_precomputations, 0L)
    expect_false(changed@functs_called[["run_qc"]])
    expect_error(get_result(changed, "run_discovery_analysis"), "No completed current result")
    expect_warning(saved <- get_result(changed, "run_discovery_analysis", historical = TRUE), "historical")
    attr(saved, "sceptre_result_record") <- NULL
    expect_identical(saved, historical)
    expect_true(result_record_is_current(object, "run_discovery_analysis"))
    expect_error(set_response_matrix(object, counts[, -1L]), "aligned cells")
})

test_that("response cache identity checks actual ordered fitting inputs at lookup", {
    object <- component_workflow()
    object <- suppressMessages(run_discovery_analysis(object, print_progress = FALSE))
    pairs <- object@preparation@discovery_pairs_with_info
    pairs <- pairs[pairs$pass_qc, ]
    context <- prepare_analysis_context(object, "discovery_analysis", pairs)
    expect_gt(length(context$response_fits), 0L)
    expect_identical(context$fit_identity, object@fits@identity)
    bad <- object
    bad@preparation@covariate_matrix[1, 1] <- 9
    expect_length(prepare_analysis_context(bad, "discovery_analysis", pairs)$response_fits, 0L)
    bad <- object
    bad@preparation@cells_in_use <- rev(bad@preparation@cells_in_use)
    expect_length(prepare_analysis_context(bad, "discovery_analysis", pairs)$response_fits, 0L)
    bad <- object
    bad@data@identity <- "different-count-content"
    expect_length(prepare_analysis_context(bad, "discovery_analysis", pairs)$response_fits, 0L)
})

test_that("assignment and eligibility edits invalidate tests independently of response fits", {
    for (name in c("high_complement_crt", "low_complement_perm")) {
        object <- component_workflow(name)
        object <- suppressMessages(run_discovery_analysis(object, print_progress = FALSE))
        cached <- object@fits@response_precomputations
        changed <- suppressMessages(assign_grnas(object, method = "thresholding",
            threshold = 41, print_progress = FALSE))
        expect_false(result_record_is_current(changed, "run_discovery_analysis"))
        expect_true(changed@preparation@stages$assignment$complete)
        expect_false(changed@preparation@stages$cell_qc$complete)
        if (object@data@low_moi) expect_length(changed@fits@response_precomputations, 0L)
        else expect_identical(changed@fits@response_precomputations, cached)
        changed <- suppressMessages(run_qc(object, response_n_umis_range = c(0, 1),
            response_n_nonzero_range = c(0, 1), n_nonzero_trt_thresh = 8L,
            n_nonzero_cntrl_thresh = 8L,
            additional_cells_to_remove = setdiff(seq_len(ncol(get_response_matrix(object))),
                object@preparation@cells_in_use)))
        expect_identical(changed@preparation@cells_in_use, object@preparation@cells_in_use)
        expect_identical(changed@fits@response_precomputations, cached)
        expect_false(result_record_is_current(changed, "run_discovery_analysis"))
        expect_true(changed@preparation@stages$pair_qc$complete)
    }
})


test_that("replacement refreshes imported guide assignment and response QC summaries", {
    object <- component_workflow("low_complement_perm")
    guides <- get_grna_matrix(object)
    guides@x[] <- 0
    changed <- set_grna_matrix(object, guides)
    expect_true(all(changed@data@covariate_data_frame$grna_n_umis == 0))
    expect_true(all(!is.finite(changed@data@import_grna_assignment_info$max_grna_frac_umis)))
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    changed <- suppressMessages(set_analysis_parameters(changed,
        discovery_pairs = inputs$discovery_pairs, formula_object = ~ batch + z,
        control_group = "complement", resampling_mechanism = "permutations"))
    changed <- suppressMessages(assign_grnas(changed, method = "maximum", print_progress = FALSE))
    expect_setequal(changed@preparation@cells_w_zero_or_twoplus_grnas,
        seq_len(ncol(guides)))
    responses <- get_response_matrix(object)
    responses@x[] <- responses@x + 1
    changed <- set_response_matrix(object, as(responses, "CsparseMatrix"))
    expect_s4_class(get_response_matrix(changed), "dgRMatrix")
    expect_equal(changed@data@covariate_data_frame$response_n_umis,
        as.numeric(Matrix::colSums(responses)))
    expect_false(identical(changed@data@covariate_data_frame$response_n_umis,
        object@data@covariate_data_frame$response_n_umis))
})


test_that("preparation cannot publish completion before computation succeeds", {
    object <- component_workflow()
    for (stage in c("set_analysis_parameters", "assign_grnas", "run_qc")) {
        before <- serialize(object, NULL)
        prepared <- perform_status_check_and_update(object, stage)
        expect_false(prepared@functs_called[[stage]])
        expect_identical(serialize(object, NULL), before)
    }
    prepared <- suppressMessages(run_qc_pt_1(object,
        response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0, 1)))
    expect_true(prepared@preparation@stages$cell_qc$complete)
    expect_false(prepared@functs_called[["run_qc"]])
    expect_false(prepared@preparation@stages$pair_qc$complete)
    completed <- run_qc_pt_2(prepared)
    expect_true(completed@functs_called[["run_qc"]])
    expect_true(completed@preparation@stages$pair_qc$complete)
})
