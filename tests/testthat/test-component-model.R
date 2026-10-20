component_workflow <- function(name = "high_complement_crt") {
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    case <- readRDS(test_path("fixtures", "refactor", "cases.rds"))[[name]]
    suppressMessages(refactor_prepare(inputs, case))
}

test_that("schema3 coordinator and choices contain no execution or duplicate authorities", {
    object <- component_workflow()
    expect_identical(slotNames(object), c("data", "spec", "preparation", "fits", "results", "schema_version"))
    expect_false(any(c("functs_called", "last_function_called", "nuclear", "nf_pipeline", "elements_to_analyze") %in% slotNames(object)))
    expect_false(any(c("requested", "B1", "B2", "B3", "negative_control_pairs", "calibration_group_size", "n_calibration_pairs") %in% slotNames(object@spec)))
    expect_false(any(c("covariate_names", "mitochondrial_gene", "integer_id") %in% slotNames(object@data)))
    expect_false(any(c("n_discovery_pairs", "n_ok_discovery_pairs", "mean_cells_per_grna", "grna_target_data_frame") %in% slotNames(object@preparation)))
    expect_identical(slotNames(object@results), c("records", "import_metadata"))
    expect_true(all(vapply(object@preparation@stages, function(x) identical(names(x), "identity"), logical(1))))
    expect_true(all(workflow_status(object)[1:4]))
    expect_true(validObject(object))
})

test_that("canonical requested families and metadata are not expanded or edited by preparation", {
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    for (name in c("low_bonferroni", "low_singleton")) {
        object <- component_workflow(name)
        expect_identical(object@spec@discovery_pairs$grna_target, inputs$discovery_pairs$grna_target)
        expect_false("grna_group" %in% names(object@spec@discovery_pairs))
        expect_false("grna_group" %in% names(object@data@grna_target_data_frame))
        expect_true("grna_group" %in% names(object@preparation@discovery_pairs))
        expect_equal(pair_count(object, "discovery"), nrow(inputs$discovery_pairs))
        expect_identical(guide_group_metadata(object)$grna_group, object@preparation@guide_groups)
    }
})

test_that("stages and reports agree with the single authoritative preparation record", {
    object <- component_workflow()
    object <- suppressMessages(run_discovery_analysis(object, print_progress = FALSE))
    expect_true(result_record_is_current(object, "run_discovery_analysis"))
    incomplete <- object
    incomplete@preparation@stages$assignment <- NULL
    expect_false(stage_complete(incomplete, "assignment"))
    expect_false(stage_complete(incomplete, "run_qc"))
    expect_false(workflow_status(incomplete)[["run_discovery_analysis"]])
    expect_error(get_result(incomplete, "run_discovery_analysis"), "No completed current")
    partial <- suppressMessages(run_qc_pt_1(object, response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1)))
    expect_true(stage_complete(partial, "cell_qc"))
    expect_false(stage_complete(partial, "pair_qc"))
    expect_false(workflow_status(partial)[["run_qc"]])
    complete <- run_qc_pt_2(partial)
    expect_true(stage_complete(complete, "run_qc"))
    expect_false(result_record_is_current(complete, "run_discovery_analysis"))
})

test_that("count replacement refreshes summaries and invalidates all derived stages", {
    object <- component_workflow()
    object <- suppressMessages(run_discovery_analysis(object, print_progress = FALSE))
    historical <- get_result(object, "run_discovery_analysis")
    counts <- get_response_matrix(object)
    counts@x <- counts@x + 1
    changed <- set_response_matrix(object, as(counts, "CsparseMatrix"))
    expect_s4_class(get_response_matrix(changed), "dgRMatrix")
    expect_equal(changed@data@covariate_data_frame$response_n_umis, as.numeric(Matrix::colSums(counts)))
    expect_false(identical(changed@data@identity, object@data@identity))
    expect_length(changed@fits@response_fits, 0L)
    expect_identical(changed@preparation@stages, list())
    expect_false(workflow_status(changed)[["set_analysis_parameters"]])
    expect_error(get_result(changed, "run_discovery_analysis"), "No completed current")
    expect_warning(saved <- get_result(changed, "run_discovery_analysis", historical = TRUE), "historical")
    attr(saved, "sceptre_result_record") <- NULL
    expect_identical(saved, historical)
    expect_true(result_record_is_current(object, "run_discovery_analysis"))
})

test_that("response cache rejects malformed coefficients, dispersion and input scopes", {
    object <- component_workflow()
    object <- suppressMessages(run_discovery_analysis(object, print_progress = FALSE))
    cache <- object@fits
    expect_s4_class(cache, "FitCache")
    expect_gt(length(cache@response_fits), 0L)
    expect_error(new("FitCache", response_fits = cache@response_fits, response_identity = list()), "ordered fitting identity")
    malformed <- cache@response_fits
    malformed[[1]]$theta <- NA_real_
    expect_error(new("FitCache", response_fits = malformed, response_identity = cache@response_identity), "dispersion")
    malformed <- cache@response_fits
    names(malformed[[1]]$fitted_coefs) <- rev(names(malformed[[1]]$fitted_coefs))
    expect_error(new("FitCache", response_fits = malformed, response_identity = cache@response_identity), "coefficients")
    malformed <- cache@response_fits
    names(malformed)[1] <- "nonexistent-response"
    expect_error(new("sceptre_object", data = object@data, spec = object@spec,
        preparation = object@preparation, fits = new("FitCache", response_fits = malformed,
            response_identity = cache@response_identity)), "identifiers must exist")
    pairs <- object@preparation@discovery_pairs_with_info
    pairs <- dplyr::filter(pairs, pass_qc)
    expect_gt(length(prepare_analysis_context(object, "discovery_analysis", pairs)$response_fits), 0L)
    changed <- object
    changed@preparation@design_matrix[1, 1] <- 9
    expect_length(prepare_analysis_context(changed, "discovery_analysis", pairs)$response_fits, 0L)
    entries <- cache@response_fits
    attr(entries, "fit_identity") <- cache@response_identity
    expect_error(commit_response_fits(changed, entries), "actual fitting inputs")
})

test_that("assignment and pair eligibility invalidation preserve compatible response fits", {
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    for (name in c("high_complement_crt", "low_complement_perm")) {
        object <- component_workflow(name)
        object <- suppressMessages(run_discovery_analysis(object, print_progress = FALSE))
        cached <- object@fits@response_fits
        changed <- suppressMessages(assign_grnas(object, method = "thresholding", threshold = 41, print_progress = FALSE))
        expect_true(stage_complete(changed, "assignment"))
        expect_false(stage_complete(changed, "cell_qc"))
        expect_false(result_record_is_current(changed, "run_discovery_analysis"))
        if (object@data@low_moi) expect_length(changed@fits@response_fits, 0L)
        else expect_identical(changed@fits@response_fits, cached)
        changed <- suppressMessages(run_qc(object, response_n_umis_range = c(0, 1),
            response_n_nonzero_range = c(0, 1), n_nonzero_trt_thresh = 8L, n_nonzero_cntrl_thresh = 8L,
            additional_cells_to_remove = inputs$remove_cells))
        expect_identical(changed@preparation@analysis_cells, object@preparation@analysis_cells)
        expect_identical(changed@fits@response_fits, cached)
        expect_true(stage_complete(changed, "run_qc"))
        expect_false(result_record_is_current(changed, "run_discovery_analysis"))
    }
})
