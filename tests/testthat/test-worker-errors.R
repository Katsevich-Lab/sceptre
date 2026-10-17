test_that("worker errors retain their cause before result collection", {
    run <- refactor_private("execute_analysis_blocks")
    worker <- function(block, proc_id) {
        if (block == "gene_X") stop("Cannot read response gene_X")
        list(result = block)
    }
    expect_error(run(list("gene_X", "gene_Y"), worker, FALSE),
                 "Cannot read response gene_X")
    expect_identical(run(list("gene_Y", "gene_Z"), worker, FALSE),
                     list(list(result = "gene_Y"), list(result = "gene_Z")))
    skip_on_os("windows")
    expect_warning(
        expect_error(run(list("gene_X", "gene_Y"), worker, TRUE),
                     "Analysis worker 1 failed: Cannot read response gene_X"),
        "scheduled core 1 encountered error")
    expect_identical(run(list("gene_Y", "gene_Z"), worker, TRUE),
                     list(list(result = "gene_Y"), list(result = "gene_Z")))
})

test_that("successful empty worker results remain valid", {
    run <- refactor_private("execute_analysis_blocks")
    worker <- function(block, proc_id) NULL
    expect_identical(run(list(1L), worker, FALSE), list(NULL))
    skip_on_os("windows")
    expect_identical(run(list(1L, 2L), worker, TRUE), list(NULL, NULL))
})

test_that("a child exiting without replying cannot silently drop a block", {
    skip_on_os("windows")
    run <- refactor_private("execute_analysis_blocks")
    worker <- function(block, proc_id) {
        # Terminate only this test-created child, without inherited R unload hooks.
        if (proc_id == 1L) tools::pskill(Sys.getpid(), signal = 9L)
        list(result = block)
    }
    expect_warning(
        expect_error(run(list(1L, 2L), worker, TRUE), "worker 1 returned no result"),
        "did not deliver")
})

test_that("CRT assembly retains failed-QC rows when no requests are eligible", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$high_complement_crt
    object <- refactor_capture(refactor_prepare(inputs, case))$value
    pairs <- object@discovery_pairs_with_info
    pairs$pass_qc[] <- FALSE
    run <- refactor_private("analyze_association_pairs")
    serial <- refactor_capture(run(object, pairs, "discovery_analysis", 1L, FALSE, FALSE, 2L, tempdir()))
    expect_identical(serial$warnings, character())
    expect_equal(nrow(serial$value$result), nrow(pairs))
    expect_true(all(!serial$value$result$pass_qc))
    expect_true(all(is.na(serial$value$result$p_value)))
    skip_on_os("windows")
    parallel <- refactor_capture(run(object, pairs, "discovery_analysis", 1L, FALSE, TRUE, 2L, tempdir()))
    expect_identical(parallel$warnings, character())
    expect_equal(parallel$value, serial$value)
})
