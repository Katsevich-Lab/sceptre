# Direct matrix oracle intentionally does not call the CellView/index helpers.
wilcoxon_oracle <- function(inputs, object, side) {
    cells <- setdiff(seq_len(ncol(inputs$response_matrix)), inputs$remove_cells)
    metadata <- inputs$grna_target_data_frame
    guide <- inputs$low_grna_matrix >= 40
    nt <- intersect(cells, which(colSums(guide[metadata$grna_target == "non-targeting", , drop = FALSE]) > 0))
    depth <- colSums(inputs$response_matrix)
    requests <- object@preparation@discovery_pairs_with_info
    requests$p_value <- NA_real_
    for (i in which(requests$pass_qc)) {
        id <- requests$grna_group[i]
        rows <- if (object@spec@grna_integration_strategy == "union") {
            which(metadata$grna_target == id)
        } else match(id, metadata$grna_id)
        trt <- intersect(cells, which(colSums(guide[rows, , drop = FALSE]) > 0))
        y <- 10000 * inputs$response_matrix[requests$response_id[i], ] / depth
        requests$p_value[i] <- stats::wilcox.test(log1p(y[trt]), log1p(y[nt]),
            alternative = side, exact = FALSE, correct = TRUE)$p.value
    }
    requests
}

test_that("Wilcoxon uses aligned full-library ranks and complete-family correction", {
    inputs <- refactor_make_inputs()
    for (integration in c("union", "singleton", "bonferroni")) {
        for (side in c("left", "both", "right")) {
            case <- refactor_case_definitions()$low_nt_perm
            case$integration <- integration
            case$side <- side
            object <- suppressMessages(refactor_prepare(inputs, case))
            before <- serialize(object, NULL)
            set.seed(378)
            rng <- .Random.seed
            alternative <- run_wilcoxon(object)
            expect_identical(.Random.seed, rng)
            expect_identical(serialize(object, NULL), before)
            expect_identical(alternative@fits, object@fits)
            expected <- wilcoxon_oracle(inputs, object,
                c(left = "less", both = "two.sided", right = "greater")[[side]])
            if (integration == "bonferroni") {
                expected$grna_target <- inputs$grna_target_data_frame$grna_target[
                    match(expected$grna_group, inputs$grna_target_data_frame$grna_id)]
                expected <- do.call(rbind, lapply(split(expected, paste(expected$response_id, expected$grna_target)), function(x) {
                    data.frame(response_id = x$response_id[1], grna_target = x$grna_target[1],
                        p_value = if (any(x$pass_qc)) min(1, min(x$p_value, na.rm = TRUE) * sum(x$pass_qc)) else NA_real_,
                        pass_qc = any(x$pass_qc))
                }))
            } else {
                names(expected)[names(expected) == "grna_group"] <- if (integration == "union") "grna_target" else "grna_id"
            }
            actual <- get_result(alternative, "run_discovery_analysis")
            keys <- c("response_id", if (integration == "singleton") "grna_id" else "grna_target")
            key <- function(x) do.call(paste, as.data.frame(x)[keys])
            actual <- actual[match(key(expected), key(actual)), ]
            expect_equal(actual$p_value, expected$p_value, tolerance = 1e-14)
            expect_identical(actual$pass_qc, expected$pass_qc)
            expect_identical(actual$significant, stats::p.adjust(actual$p_value,
                method = object@spec@multiple_testing_method) < object@spec@multiple_testing_alpha)
            expect_true(all(is.na(actual$p_value[!actual$pass_qc])))
            expect_true(all(actual$inference_method == "wilcoxon_unadjusted"))
            expect_true(all(actual$effect_estimator == "mean_CP10K_ratio"))
            record <- alternative@results@records$run_discovery_analysis
            expect_false(record$provenance$covariate_adjustment)
            expect_identical(record$provenance$response_fit_method, "not_used")
            expect_identical(record$provenance$plan@B1, 0L)
            expect_identical(record$provenance$plan@B2, 0L)
            expect_identical(record$provenance$plan@B3, 0L)
            expect_true(all(is.na(actual$se_fold_change)))
            path <- withr::local_tempfile()
            saveRDS(alternative, path)
            reopened <- readRDS(path)
            expect_equal(get_result(reopened, "run_discovery_analysis"), get_result(alternative, "run_discovery_analysis"))
        }
    }
})

test_that("unsupported Wilcoxon modes reject before execution", {
    expect_error(run_wilcoxon(NULL, execution = "distributed"), "direct serial")
    expect_error(run_wilcoxon(NULL, parallel = TRUE), "direct serial")
    inputs <- refactor_make_inputs()
    object <- suppressMessages(refactor_prepare(inputs, refactor_case_definitions()$low_complement_perm))
    expect_error(run_wilcoxon(object), "low MOI and NT")
    object <- suppressMessages(refactor_prepare(inputs, refactor_case_definitions()$high_complement_crt))
    expect_error(run_wilcoxon(object), "low MOI and NT")
    object <- suppressMessages(refactor_prepare(inputs, refactor_case_definitions()$low_nt_perm))
    alternative <- run_wilcoxon(object)
    expect_error(prepare_execution(alternative, "run_discovery_analysis"), "Wilcoxon is unsupported")
    expect_error(run_wilcoxon(object, n_calibration_pairs = 2), "Calibration sizes")
})

test_that("all-tied ranks, null calibration and method transitions are explicit", {
    for (side in c("less", "two.sided", "greater")) {
        expect_identical(sceptre:::wilcoxon_probability(rep(3, 12), rep(3, 15), side), 1)
    }
    inputs <- refactor_make_inputs()
    object <- suppressMessages(refactor_prepare(inputs, refactor_case_definitions()$low_nt_perm))
    alternative <- suppressMessages(run_wilcoxon(object, "calibration_check", n_calibration_pairs = 4L, calibration_group_size = 1L))
    expect_equal(nrow(get_result(alternative, "run_calibration_check")), 4L)
    alternative <- run_wilcoxon(alternative, "power_check")
    alternative <- run_wilcoxon(alternative)
    expect_true(sceptre:::result_record_is_current(alternative, "run_calibration_check"))
    expect_output(print(alternative), "experimental unadjusted Wilcoxon")
    path <- withr::local_tempfile()
    sceptre:::write_analysis_summary(alternative, path)
    expect_match(paste(readLines(path, warn = FALSE), collapse = "\n"), "no covariate adjustment")
    # A public SCEPTRE entry must invalidate Wilcoxon before consulting calibration.
    # Run power on this tiny fixture; no changed resampling policy is used.
    standard <- suppressMessages(run_power_check(alternative, print_progress = FALSE))
    expect_error(get_result(standard, "run_calibration_check"), "No completed current")
    expect_error(get_result(standard, "run_discovery_analysis"), "No completed current")
    expect_warning(historical <- get_result(standard, "run_discovery_analysis", historical = TRUE), "historical")
    expect_identical(historical$inference_method,
        get_result(alternative, "run_discovery_analysis")$inference_method)
    again <- run_wilcoxon(standard)
    expect_error(get_result(again, "run_power_check"), "No completed current")
    expect_equal(get_result(again, "run_discovery_analysis")$p_value,
        get_result(alternative, "run_discovery_analysis")$p_value)
})

test_that("eligible all-tied families return one and failed-QC families remain missing", {
    inputs <- refactor_make_inputs()
    inputs$response_matrix[,] <- 3
    object <- suppressMessages(refactor_prepare(inputs, refactor_case_definitions()$low_nt_perm))
    result <- get_result(run_wilcoxon(object), "run_discovery_analysis")
    expect_true(all(result$p_value == 1))
    expect_true(all(result$fold_change == 1))
    expect_false(any(result$significant))
    inputs$response_matrix[,] <- 0
    # Preserve positive library depth, but request only zero genes.
    inputs$response_matrix[8L, ] <- 20
    object <- suppressMessages(refactor_prepare(inputs, refactor_case_definitions()$low_nt_perm))
    result <- get_result(run_wilcoxon(object), "run_discovery_analysis")
    expect_equal(nrow(result), nrow(inputs$discovery_pairs))
    expect_true(all(is.na(result$p_value)))
    expect_true(all(result$test_status == "failed_qc"))
})

test_that("alternative plots label the actual method and descriptive effect", {
    inputs <- refactor_make_inputs()
    object <- suppressMessages(refactor_prepare(inputs, refactor_case_definitions()$low_nt_perm))
    alternative <- run_wilcoxon(object)
    before <- serialize(alternative, NULL)
    expect_no_warning(p <- plot_run_discovery_analysis(alternative))
    expect_identical(serialize(alternative, NULL), before)
    expect_match(p$labels$caption, "Experimental unadjusted Wilcoxon")
    plots <- plot_run_discovery_analysis(alternative, return_indiv_plots = TRUE)
    expect_true(all(vapply(plots, function(p) grepl("mean CP10K ratio", p$labels$caption), logical(1))))
})

test_that("Wilcoxon calibration never computes an unused resampling budget", {
    object <- suppressMessages(refactor_prepare(refactor_make_inputs(), refactor_case_definitions()$low_nt_perm))
    testthat::local_mocked_bindings(run_resampling_budget = function(...) stop("unused budget"), .package = "sceptre")
    result <- suppressMessages(run_wilcoxon(object, "calibration_check", n_calibration_pairs = 4L, calibration_group_size = 1L))
    expect_equal(nrow(get_result(result, "run_calibration_check")), 4L)
})

test_that("public method accessor reports current families and rejects malformed methods", {
    object <- suppressMessages(refactor_prepare(refactor_make_inputs(), refactor_case_definitions()$low_nt_perm))
    expect_identical(get_result_methods(object), stats::setNames(character(), character()))
    alternative <- run_wilcoxon(object)
    expect_identical(get_result_methods(alternative), c(run_discovery_analysis = "wilcoxon_unadjusted"))
    old <- alternative
    old@results@records$run_discovery_analysis$invalidated_by <- "stale"
    expect_identical(get_result_methods(old), stats::setNames(character(), character()))
    bad <- alternative
    bad@results@records$run_discovery_analysis$provenance$inference_method <- NA_character_
    expect_error(get_result_methods(bad), "Invalid inference method")
    old@results@records$run_discovery_analysis$invalidated_by <- NULL
    old@results@records$run_discovery_analysis$provenance$inference_method <- NULL
    expect_identical(get_result_methods(old), c(run_discovery_analysis = "sceptre"))
    old@spec@side_code <- -1L
    expect_identical(get_result_methods(old), stats::setNames(character(), character()))
})

test_that("alternative display floor preserves scientific probabilities", {
    raw <- data.frame(p_value = c(0, NA_real_, 0.04), inference_method = "wilcoxon_unadjusted")
    display <- sceptre:::alternative_plot_probabilities(raw)
    expect_identical(display$p_value, c(1e-300, NA_real_, 0.04))
    expect_identical(raw$p_value, c(0, NA_real_, 0.04))
    standard <- raw
    standard$inference_method <- NULL
    expect_identical(sceptre:::alternative_plot_probabilities(standard), standard)
})

test_that("a pre-existing SCEPTRE manifest cannot bypass the Wilcoxon guard", {
    object <- suppressMessages(refactor_prepare(refactor_make_inputs(), refactor_case_definitions()$low_nt_perm))
    manifest <- prepare_execution(object, "run_discovery_analysis")
    alternative <- run_wilcoxon(object)
    testthat::local_mocked_bindings(execute_unit_analysis = function(...) stop("worker must not run"), .package = "sceptre")
    expect_error(execute_work_unit(alternative, manifest, manifest$units[[1L]]$id), "Wilcoxon is unsupported")
    expect_error(collect_execution(alternative, manifest, list()), "Wilcoxon is unsupported")
})

test_that("standard tibble plots do not inspect an absent alternative column", {
    standard <- dplyr::tibble(p_value = 0.04)
    expect_no_warning(value <- sceptre:::alternative_volcano_threshold(0.04, standard))
    expect_identical(value, 0.04)
})

test_that("Wilcoxon drops inherited fitter-native usage provenance", {
    object <- suppressMessages(refactor_prepare(refactor_make_inputs(), refactor_case_definitions()$low_nt_perm))
    inherited <- sceptre:::completed_result_provenance
    testthat::local_mocked_bindings(completed_result_provenance = function(...) {
        provenance <- inherited(...)
        provenance$perturbation_fit_native <- "unused-selected-fitter-native"
        provenance
    }, .package = "sceptre")
    alternative <- run_wilcoxon(object)
    provenance <- alternative@results@records$run_discovery_analysis$provenance
    expect_false("perturbation_fit_native" %in% names(provenance))
    expect_identical(provenance$response_fit_method, "not_used")
    expect_identical(provenance$perturbation_fit_method, "not_used")
    expect_identical(get_result_methods(alternative), c(run_discovery_analysis = "wilcoxon_unadjusted"))
    expect_match(plot_run_discovery_analysis(alternative)$labels$caption,
        "Experimental unadjusted Wilcoxon")
})

test_that("selected fast fitter remains unused by combined-API Wilcoxon", {
    skip_if_not("set_perturbation_fit_method" %in% getNamespaceExports("sceptre"),
        "W07 public fitter selector is not present in this standalone branch")
    object <- suppressMessages(refactor_prepare(refactor_make_inputs(), refactor_case_definitions()$low_nt_perm))
    select_fitter <- getExportedValue("sceptre", "set_perturbation_fit_method")
    get_fitter <- getExportedValue("sceptre", "get_perturbation_fit_method")
    selected <- select_fitter(object, "fast_logistic")
    testthat::local_mocked_bindings(
        fit_perturbation_model = function(...) stop("unused fitter must not execute"),
        fit_response = function(...) stop("unused response fit must not execute"),
        .package = "sceptre")
    expected <- get_result(run_wilcoxon(object), "run_discovery_analysis")
    alternative <- run_wilcoxon(selected)
    provenance <- alternative@results@records$run_discovery_analysis$provenance
    expect_false("perturbation_fit_native" %in% names(provenance))
    expect_identical(provenance$response_fit_method, "not_used")
    expect_identical(provenance$perturbation_fit_method, "not_used")
    expect_identical(get_result_methods(alternative), c(run_discovery_analysis = "wilcoxon_unadjusted"))
    expect_identical(get_fitter(selected), "fast_logistic")
    expect_equal(get_result(alternative, "run_discovery_analysis")$p_value,
        expected$p_value)
})
