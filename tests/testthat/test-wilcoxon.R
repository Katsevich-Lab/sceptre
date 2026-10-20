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
    inputs$response_matrix[1L, ] <- inputs$response_matrix[5L, ]
    object <- suppressMessages(refactor_prepare(inputs, refactor_case_definitions()$low_nt_perm))
    alternative <- run_wilcoxon(object)
    p <- plot_run_discovery_analysis(alternative)
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
