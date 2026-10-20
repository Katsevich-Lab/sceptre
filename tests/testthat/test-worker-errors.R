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

test_that("CRT dispatch preserves empty internal request blocks", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$high_complement_crt
    object <- refactor_capture(refactor_prepare(inputs, case))$value
    pairs <- object@preparation@discovery_pairs_with_info[FALSE, , drop = FALSE]
    context <- refactor_private("prepare_analysis_context")(object, "discovery_analysis", pairs)
    execution <- list(output_amount = 1L, print_progress = FALSE, parallel = FALSE,
                      n_processors = 2L, log_dir = tempdir())
    run <- refactor_private("run_crt_analysis")
    serial <- refactor_capture(run(context, execution))
    expect_identical(serial$warnings, character())
    expect_equal(nrow(serial$value$result), 0L)
    expect_identical(serial$value$response_precomputations, context$response_fits)
    skip_on_os("windows")
    execution$parallel <- TRUE
    parallel <- refactor_capture(run(context, execution))
    expect_identical(parallel$warnings, character())
    expect_equal(parallel$value, serial$value)
})

test_that("parallel mixture assignment surfaces worker cause before flattening", {
    skip_on_os("windows")
    ids <- paste0("guide-", seq_len(4L))
    counts <- matrix(0, 4L, 5L, dimnames = list(ids, NULL))
    hyperparameters <- list(n_em_rep = 2L, pi_guess_range = c(1e-5, 0.1),
        g_pert_guess_range = c(10, 50), formula_object = ~ 1,
        n_nonzero_cells_cutoff = 6L, backup_threshold = 20,
        probability_threshold = 0.8)
    testthat::local_mocked_bindings(load_row = function(...) {
        stop("guide backing store is unreadable")
    }, .package = "sceptre")
    run <- refactor_private("assign_grnas_to_cells_mixture")
    expect_warning(expect_error(suppressMessages(run(counts,
        data.frame(z = seq_len(5L)), hyperparameters, FALSE, TRUE, 2L,
        tempdir(), ids)), "Analysis worker 1 failed: guide backing store is unreadable"),
        "encountered errors")
})
