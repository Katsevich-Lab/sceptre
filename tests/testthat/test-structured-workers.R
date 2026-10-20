test_that("serial and fork failures preserve structured causes and execution context", {
    run <- refactor_private("execute_analysis_blocks")
    cause <- structure(list(message = "fixture backing failure", call = NULL,
        store = "response"), class = c("fixture_store_error", "error", "condition"))
    worker <- function(block, proc_id) stop(cause)
    observe <- function(parallel) suppressWarnings(tryCatch(
        run(list("r1", "r2"), worker, parallel, "discovery_analysis", "response_fits"),
        sceptre_worker_error = identity))
    check <- function(error) {
        expect_s3_class(error, "sceptre_worker_error")
        expect_identical(error$block, 1L)
        expect_identical(error$analysis, "discovery_analysis")
        expect_identical(error$phase, "response_fits")
        expect_identical(error$reason, "error")
        expect_identical(error$parent, cause)
        expect_identical(error$cause, cause)
        expect_match(conditionMessage(error), "Analysis worker 1 failed: fixture backing failure", fixed = TRUE)
    }
    check(observe(FALSE))
    skip_on_os("windows")
    check(observe(TRUE))
})

test_that("worker interrupts escape unchanged in serial and fork execution", {
    run <- refactor_private("execute_analysis_blocks")
    interruption <- structure(list(message = "fixture interrupted", call = NULL),
        class = c("fixture_interrupt", "interrupt", "condition"))
    worker <- function(block, proc_id) stop(interruption)
    observe <- function(parallel) tryCatch(run(list(1L, 2L), worker, parallel), interrupt = identity)
    expect_identical(observe(FALSE), interruption)
    skip_on_os("windows")
    expect_identical(observe(TRUE), interruption)
})

test_that("lost child replies have a machine-readable reason", {
    skip_on_os("windows")
    run <- refactor_private("execute_analysis_blocks")
    failure <- suppressWarnings(tryCatch(run(list(1L, 2L), function(block, proc_id) {
        if (proc_id == 1L) tools::pskill(Sys.getpid(), signal = 9L)
        NULL
    }, TRUE, phase = "guide_tests"), sceptre_worker_error = identity))
    expect_identical(failure$reason, "lost_reply")
    expect_identical(failure$block, 1L)
    expect_identical(failure$phase, "guide_tests")
    expect_null(failure$parent)
})
