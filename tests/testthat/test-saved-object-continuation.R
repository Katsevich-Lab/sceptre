# Physical fixtures were serialized by their historical package definitions,
# rather than assembled under today's S4 classes. Keep their original bytes.
w02_saved_paths <- function() test_path("fixtures", "saved-objects", c(
    "sceptre-0.99.0-21f9ea0.rds", "sceptre-schema2-4aad80b.rds",
    "sceptre-schema3-list-plans-2a0938e.rds"))

test_that("physical historical layouts inspect, continue, save and reopen", {
    paths <- w02_saved_paths()
    hashes <- tools::md5sum(paths)
    for (path in paths) {
        old <- readRDS(path)
        before <- serialize(old, NULL)
        # Accessors and printing must inspect the physical old representation
        # without changing its attributes or blessing historical results.
        expect_no_error(capture.output(show(old)))
        counts <- get_response_matrix(old)
        guides <- get_grna_matrix(old)
        object <- sceptre:::migrate_sceptre_object(old)
        if (any(object@results@import_metadata$migrated_from_schema %in% c(1L, 2L))) {
            expect_false(sceptre:::result_record_is_current(object, "run_discovery_analysis"))
        }
        set.seed(2371L)
        continued <- suppressMessages(run_discovery_analysis(object, print_progress = FALSE))
        expect_true(sceptre:::result_record_is_current(continued, "run_discovery_analysis"))
        expect_identical(get_response_matrix(continued), counts)
        expect_identical(get_grna_matrix(continued), guides)
        expect_identical(serialize(old, NULL), before)
        saved <- tempfile(fileext = ".rds")
        saveRDS(continued, saved)
        reopened <- readRDS(saved)
        expect_true(validObject(reopened))
        expect_equal(get_result(reopened, "run_discovery_analysis"),
                     get_result(continued, "run_discovery_analysis"))
        rerun <- suppressMessages(run_discovery_analysis(reopened, print_progress = FALSE))
        expect_equal(get_result(rerun, "run_discovery_analysis"),
                     get_result(continued, "run_discovery_analysis"))
        expect_true(sceptre:::result_record_is_current(rerun, "run_discovery_analysis"))
        unlink(saved)
    }
    expect_identical(tools::md5sum(paths), hashes)
})

test_that("failed or interrupted historical continuation leaves original recoverable", {
    old <- readRDS(w02_saved_paths()[[3L]])
    before <- serialize(old, NULL)
    original <- sceptre:::analyze_association_pairs
    interruption <- structure(list(message = "controlled continuation interrupt", call = NULL),
        class = c("w02_interrupt", "interrupt", "condition"))
    for (kind in c("error", "interrupt")) {
        testthat::local_mocked_bindings(analyze_association_pairs = function(...) {
            if (kind == "interrupt") stop(interruption)
            stop("controlled continuation failure")
        }, .package = "sceptre")
        failure <- tryCatch(suppressMessages(run_discovery_analysis(old, print_progress = FALSE)),
            error = identity, interrupt = identity)
        expect_s3_class(failure, if (kind == "error") "error" else "w02_interrupt")
        expect_identical(serialize(old, NULL), before)
    }
    testthat::local_mocked_bindings(analyze_association_pairs = original, .package = "sceptre")
    recovered <- suppressMessages(run_discovery_analysis(old, print_progress = FALSE))
    expect_true(sceptre:::result_record_is_current(recovered, "run_discovery_analysis"))
    expect_identical(serialize(old, NULL), before)
})
