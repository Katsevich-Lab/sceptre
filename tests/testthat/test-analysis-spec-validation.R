test_that("correction validation preserves documented endpoints and stats aliases", {
    for (alpha in c(0.01, 0.1, 1)) {
        for (method in c(stats::p.adjust.methods, "f", "bo", "hol", "hoch")) {
            expect_silent(validate_correction_settings(alpha, method))
            spec <- new("AnalysisSpec", multiple_testing_alpha = alpha,
                multiple_testing_method = method)
            expect_true(validObject(spec))
            expect_identical(spec@multiple_testing_method, method)
        }
    }
    expect_identical(compute_empirical_budget(3L, 1, 0L, 499L, 0L), 30L)
    expect_true(validObject(new("AnalysisSpec")))
})

test_that("malformed correction values fail before default formula work", {
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    case <- readRDS(test_path("fixtures", "refactor", "cases.rds"))$low_complement_perm
    object <- refactor_import(inputs, case)
    before <- serialize(object, NULL)
    testthat::local_mocked_bindings(auto_construct_formula_object = function(...) stop("expensive formula construction reached"))
    for (alpha in list(NULL, numeric(), c(0.05, 0.1), NA_real_, NaN, Inf, -Inf,
        0, -0.1, 1.5, TRUE, "0.1", 0.1 + 0i, matrix(0.1, 1, 1))) {
        error <- tryCatch(set_analysis_parameters(object, multiple_testing_alpha = alpha), error = identity)
        expect_s3_class(error, "sceptre_argument_error")
        expect_identical(error$argument, "multiple_testing_alpha")
        expect_identical(serialize(object, NULL), before)
    }
    for (method in list(NULL, character(), c("BH", "BY"), NA_character_, "", "h", "unknown", 1, factor("BH"), matrix("BH", 1, 1))) {
        error <- tryCatch(set_analysis_parameters(object, multiple_testing_method = method), error = identity)
        expect_s3_class(error, "sceptre_argument_error")
        expect_identical(error$argument, "multiple_testing_method")
        expect_identical(serialize(object, NULL), before)
    }
})

test_that("initialized AnalysisSpec reuses scalar correction semantics", {
    original <- new("AnalysisSpec", multiple_testing_alpha = 0.1, multiple_testing_method = "BH")
    for (alpha in list(numeric(), c(0.05, 0.1), NA_real_, Inf, 0, -1, 2)) {
        spec <- original
        spec@multiple_testing_alpha <- alpha
        expect_match(validObject(spec, test = TRUE), "multiple_testing_alpha")
    }
    for (method in list(character(), c("BH", "BY"), NA_character_, "h", "unknown")) {
        spec <- original
        spec@multiple_testing_method <- method
        expect_match(validObject(spec, test = TRUE), "multiple_testing_method")
    }
    spec <- new("AnalysisSpec")
    spec@side_code <- 0L
    expect_true(any(grepl("multiple_testing_alpha", validObject(spec, test = TRUE))))
})

test_that("valid nondefault correction and alpha one survive public setup and finalization", {
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    case <- readRDS(test_path("fixtures", "refactor", "cases.rds"))$low_complement_perm
    object <- refactor_import(inputs, case)
    object <- set_analysis_parameters(object, formula_object = ~ batch + z,
        discovery_pairs = inputs$discovery_pairs,
        multiple_testing_alpha = 1, multiple_testing_method = "bo")
    expect_identical(object@spec@multiple_testing_method, "bo")
    expect_identical(object@spec@multiple_testing_alpha, 1)
    expect_true(validObject(object@spec))
    malformed <- object
    malformed@spec@multiple_testing_alpha <- c(0.1, 0.2)
    expect_true(any(grepl("multiple_testing_alpha", validObject(malformed, test = TRUE))))
    table <- data.frame(response_id = c("r1", "r2"), grna_group = c("A", "B"), p_value = c(0.01, 0.6))
    result <- process_discovery_result(table, object)
    expect_identical(result$significant, stats::p.adjust(table$p_value, method = "bonferroni") < 1)
})
