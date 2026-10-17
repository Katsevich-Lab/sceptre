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

test_that("missing worker results cannot silently drop a block", {
    run <- refactor_private("execute_analysis_blocks")
    worker <- function(block, proc_id) NULL
    expect_error(run(list(1L), worker, FALSE), "worker 1 returned no result")
    skip_on_os("windows")
    expect_error(run(list(1L, 2L), worker, TRUE), "worker 1 returned no result")
})
