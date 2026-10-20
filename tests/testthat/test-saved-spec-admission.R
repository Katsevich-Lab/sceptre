saved_spec_roundtrip <- function(object) {
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path))
    saveRDS(object, path)
    readRDS(path)
}

saved_spec_crt_object <- function() {
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    case <- readRDS(test_path("fixtures", "refactor", "cases.rds"))$high_complement_crt
    object <- suppressMessages(refactor_prepare(inputs, case))
    suppressMessages(run_discovery_analysis(object, print_progress = FALSE,
        parallel = FALSE, output_amount = 1L))
}

test_that("public consumers reject corrupt saved working settings before warm CRT work", {
    object <- saved_spec_crt_object()
    expect_gt(length(object@fits@perturbation_fits), 0L)
    testthat::local_mocked_bindings(
        prepare_analysis_context = function(...) stop("analysis context reached"),
        fit_perturbation_glm = function(...) stop("CRT fitting reached"))
    for (field in c("multiple_testing_alpha", "multiple_testing_method")) {
        values <- if (field == "multiple_testing_alpha") list(0, NA_real_, c(0.1, 0.2)) else
            list("unknown", NA_character_, c("BH", "BY"))
        for (value in values) {
            corrupt <- object
            slot(corrupt@spec, field) <- value
            restored <- saved_spec_roundtrip(corrupt)
            before <- serialize(restored, NULL)
            for (consume in list(
                function(x) run_discovery_analysis(x, print_progress = FALSE),
                function(x) run_power_check(x, print_progress = FALSE),
                function(x) get_result(x, "run_discovery_analysis"))) {
                condition <- tryCatch(consume(restored), error = identity)
                expect_s3_class(condition, "sceptre_argument_error")
                expect_identical(condition$argument, field)
                expect_identical(serialize(restored, NULL), before)
            }
        }
    }
})

test_that("saved valid aliases and alpha one are admitted before public computation", {
    object <- saved_spec_crt_object()
    testthat::local_mocked_bindings(prepare_analysis_context = function(...) stop("admitted to context"))
    for (method in c("fdr", "f", "bo", "BH")) {
        object@spec@multiple_testing_method <- method
        object@spec@multiple_testing_alpha <- 1
        restored <- saved_spec_roundtrip(object)
        expect_identical(migrate_sceptre_object(restored)@spec@multiple_testing_method, method)
        expect_error(run_discovery_analysis(restored, print_progress = FALSE), "admitted to context")
    }
    expect_silent(validate_analysis_spec(new("AnalysisSpec")))
    expect_s4_class(migrate_sceptre_object(saved_spec_roundtrip(new("sceptre_object"))), "sceptre_object")
})

test_that("saved historical records do not require whole-object validity at admission", {
    object <- saved_spec_crt_object()
    object@results@records$run_discovery_analysis$invalidated_by <- "historical test"
    # Deliberately stale preparation is compatible with historical inspection.
    object@preparation@design_matrix <- matrix(numeric(), 0L, 0L)
    restored <- saved_spec_roundtrip(object)
    expect_warning(table <- get_result(restored, "run_discovery_analysis", historical = TRUE), "historical")
    expect_s3_class(table, "data.frame")
})

test_that("legacy saved settings retain valid correction aliases at migration", {
    fixtures <- c("sceptre-0.99.0-21f9ea0.rds", "sceptre-schema2-4aad80b.rds")
    for (version in seq_along(fixtures)) {
        for (method in c("fdr", "bo")) {
            old <- readRDS(test_path("fixtures", "saved-objects", fixtures[[version]]))
            if (version == 1L) {
                attr(old, "multiple_testing_alpha") <- 1
                attr(old, "multiple_testing_method") <- method
            } else {
                spec <- attr(old, "spec")
                attr(spec, "multiple_testing_alpha") <- 1
                attr(spec, "multiple_testing_method") <- method
                attr(old, "spec") <- spec
            }
            migrated <- migrate_sceptre_object(saved_spec_roundtrip(old))
            expect_identical(migrated@spec@multiple_testing_alpha, 1)
            expect_identical(migrated@spec@multiple_testing_method, method)
            expect_s4_class(migrate_sceptre_object(saved_spec_roundtrip(migrated)), "sceptre_object")
        }
    }
})
