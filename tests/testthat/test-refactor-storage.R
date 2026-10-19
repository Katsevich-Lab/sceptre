# Small examples spell out the current count rules, including differences
# between matrix formats that the refactor must preserve.

test_that("legacy storage retains stored-zero and fractional semantics", {
    mat <- storage_test_matrix()
    expect_identical(load_row(mat, "r1"), c(0, 0, 0, 2.5, 0, 0, 0, 3.8, 0, 0, 4.2, 0))
    expect_identical(threshold_count_matrix(mat@j, mat@p, mat@x, 1L, 2.8), c(8L, 11L))
    expect_identical(threshold_count_matrix(mat@j, mat@p, mat@x, 1L, 0), c(1L, 4L, 8L, 11L))
    # A stored zero is counted by position, even though its loaded value is zero.
    counts <- sparse_pair_counts(
        mat@j, mat@p, 12L, list(c(2L, 3L)), list(4L, c(1L, 5L)),
        c(1L, 2L), c(1L, 1L), TRUE, c(2L, 4L, 5L, 8L, 11L)
    )
    expect_identical(counts$n_nonzero_trt, c(1L, 1L))
    expect_identical(counts$n_nonzero_cntrl, c(2L, 2L))
    expect_identical(counts$n_nonzero_mat, matrix(c(1L, 1L, 1L, 1L), 2L))
    expect_identical(counts$n_nonzero_tot, c(3L, 3L))
})

test_that("owned sparse group coordinates remain unchanged across calls", {
    mat <- storage_test_matrix()
    groups <- list(c(2L, 3L), c(4L, 1L, 5L))
    before <- serialize(groups, NULL)
    out <- sparse_group_counts(mat@j, mat@p, 12L, groups,
                                      c(2L, 4L, 5L, 8L, 11L))
    expect_identical(out, matrix(c(1L, 2L, 1L, 2L), 2L))
    expect_identical(serialize(groups, NULL), before)
    expect_identical(sparse_group_counts(mat@j, mat@p, 12L, groups,
        c(2L, 4L, 5L, 8L, 11L)), out)
})

test_that("legacy matrix conversion agrees with the selected ondisc copy", {
    skip_if_not_installed("ondisc")
    mat <- storage_test_matrix()
    for (input in list(mat, methods::as(mat, "CsparseMatrix"),
                       methods::as(mat, "TsparseMatrix"), as.matrix(mat))) {
        for (row_access in c(TRUE, FALSE)) {
            # Older ondisc binaries read one pointer entry past the end in
            # their conversion helper. This comparison deliberately calls that
            # independent upstream oracle; fixed binaries emit no warning.
            known <- c("subscript out of bounds (index 8 >= vector size 8)",
                       "subscript out of bounds (index 6 >= vector size 6)")
            warnings <- character()
            reference <- withCallingHandlers(
                ondisc:::set_matrix_accessibility(input, row_access),
                warning = function(w) {
                    detail <- conditionMessage(w)
                    warnings <<- c(warnings, detail)
                    if (detail %in% known) invokeRestart("muffleWarning")
                })
            expect_true(all(warnings %in% known), info = paste(warnings, collapse = "; "))
            expect_identical(set_matrix_accessibility(input, row_access), reference)
        }
    }
})

test_that("legacy backed storage coerces fractions and retains stored positions", {
    skip_if_not_installed("ondisc")
    directory <- tempfile("storage-reference-")
    dir.create(directory)
    # ondisc's legacy tilde expander misreads Windows RUNNER~1 short names.
    directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    mat <- storage_test_matrix(200L) # exceed the legacy writer's chunk size
    backed <- ondisc:::create_odm_from_r_matrix_internal(
        mat, file.path(directory, "counts.odm"), integer_id = 1234L
    )
    expect_identical(backed@integer_id, 1234L)
    expect_equal(as.numeric(backed[1L, ])[c(1L, 4L, 8L, 11L)], c(0, 2, 3, 4))
    # The existing HDF5 threshold is integer-valued (2.8 becomes 2), whereas
    # in-memory thresholding compares against the double supplied by the user.
    expect_identical(ondisc:::threshold_count_matrix_cpp(
        backed@h5_file, backed@ptr, 1L, 2.8
    )[seq_len(3L)], c(4L, 8L, 11L))
    groups <- list(c(2L, 3L), c(4L, 1L, 5L))
    out <- ondisc:::compute_n_trt_cells_matrix_ondisc(
        backed@h5_file, backed@ptr, ncol(mat), 5L, 2L, groups,
        c(2L, 4L, 5L, 8L, 11L)
    )
    expect_identical(out, matrix(c(1L, 2L, 1L, 2L), 2L))
    expect_identical(groups, list(c(1L, 2L), c(3L, 0L, 4L)))
})

test_that("storage reads preserve requested source order and repeated cells", {
    mat <- storage_test_matrix()
    expect_identical(storage_response_index(mat, "r2"), 2L)
    expect_identical(storage_read_response(mat, 1L, c(11L, 4L, 11L, 8L, 1L)),
                     c(4.2, 2.5, 4.2, 3.8, 0))
    expect_identical(storage_read_response(mat, 2L, 8L), 0)
    expect_identical(storage_read_response(mat, 2L, integer()), numeric())
    expect_identical(storage_read_response(mat, 1L), load_row(mat, "r1"))
    expect_identical(storage_threshold_row(mat, 1L, 2.8), c(8L, 11L))
})

storage_test_contract <- function(mat) {
    cells <- c(2L, 4L, 5L, 8L, 11L)
    groups <- list(A = c(2L, 3L), B = c(4L, 1L, 5L))
    nt_groups <- list(nt1 = 4L, nt2 = c(1L, 5L))
    controls <- c(4L, 1L, 5L)
    pairs <- data.frame(response_idx = c(2L, 1L, 2L, 1L),
                        group_idx = c(2L, 1L, 1L, 1L))
    original <- serialize(list(cells, groups, nt_groups, controls, pairs), NULL)
    for (control in list("complement", controls)) {
        for (nuclear in c(FALSE, TRUE)) {
            out <- expect_no_warning(storage_pair_counts(mat, pairs, groups, cells,
                control, nt_groups, nuclear))
            expect_identical(out$n_nonzero_trt, c(2L, 1L, 1L, 1L))
            expect_identical(out$n_nonzero_cntrl,
                             if (identical(control, "complement")) c(1L, 2L, 2L, 2L)
                             else rep(2L, 4L))
            if (!(nuclear && methods::is(mat, "odm"))) {
                expect_identical(out$n_nonzero_mat, matrix(rep(1L, 4L), 2L))
                expect_identical(out$n_nonzero_tot,
                                 if (identical(control, "complement")) c(3L, 3L) else c(2L, 2L))
            }
        }
        # Pair-only API must agree without constructing unused NT summaries.
        pair_only <- expect_no_warning(storage_pair_counts(mat, pairs, groups,
            cells, control))
        expect_identical(pair_only$n_nonzero_trt, out$n_nonzero_trt)
        expect_identical(pair_only$n_nonzero_cntrl, out$n_nonzero_cntrl)
        expect_identical(pair_only$n_nonzero_tot, integer())
        empty <- expect_no_warning(storage_pair_counts(mat, pairs[integer(), ],
            groups, cells, control))
        expect_identical(empty$n_nonzero_trt, integer())
        expect_identical(empty$n_nonzero_cntrl, integer())
        empty_summary <- expect_no_warning(storage_pair_counts(mat,
            pairs[integer(), ], groups, cells, control, nt_groups))
        expect_identical(empty_summary$n_nonzero_trt, integer())
        expect_identical(empty_summary$n_nonzero_cntrl, integer())
        expect_identical(empty_summary$n_nonzero_mat,
            matrix(rep(1L, 4L), 2L))
        expect_identical(empty_summary$n_nonzero_tot,
            if (identical(control, "complement")) c(3L, 3L) else c(2L, 2L))
    }
    expected <- matrix(c(1L, 2L, 1L, 2L), 2L)
    for (repeat_call in seq_len(2L)) {
        expect_identical(storage_group_counts(mat, c(2L, 1L, 2L), groups, cells),
                         expected[, c(2L, 1L, 2L), drop = FALSE])
    }
    expect_identical(storage_group_counts(mat, 1L, groups[1L], cells),
                     matrix(1L, 1L, 1L))
    expect_identical(storage_group_counts(mat, integer(), groups, cells),
                     matrix(integer(), 2L, 0L))
    expect_identical(storage_group_counts(mat, c(1L, 2L), list(), cells),
                     matrix(integer(), 0L, 2L))
    # Reordering cells kept after QC also changes group positions within that
    # list; update both so they still select the same original matrix columns.
    reordered <- storage_pair_counts(
        mat, pairs, list(c(4L, 5L), c(1L, 2L, 3L)),
        cells[c(4L, 1L, 5L, 2L, 3L)], "complement"
    )
    expect_identical(reordered$n_nonzero_trt, c(2L, 1L, 1L, 1L))
    expect_identical(reordered$n_nonzero_cntrl, c(1L, 2L, 2L, 2L))
    expect_identical(serialize(list(cells, groups, nt_groups, controls, pairs), NULL), original)
}

test_that("sparse storage counts preserve requests and caller-owned indices", {
    storage_test_contract(storage_test_matrix())
})

test_that("backed storage has the same ordered count and ownership contracts", {
    skip_if_not_installed("ondisc")
    directory <- tempfile("storage-adapter-")
    dir.create(directory)
    directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    mat <- storage_test_matrix(200L)
    set.seed(79)
    rng <- .Random.seed
    pair <- storage_create_backed_pair(mat, mat, directory, 9876L)
    expect_identical(.Random.seed, rng)
    expect_identical(pair$response_matrix@integer_id, 9876L)
    expect_identical(pair$grna_matrix@integer_id, 9876L)
    expect_identical(storage_read_response(pair$response_matrix, 1L, c(11L, 4L, 11L, 8L, 1L)),
                     c(4, 2, 4, 3, 0))
    storage_test_contract(pair$response_matrix)
    reopened <- storage_reopen_backed_pair(file.path(directory, "response.odm"),
                                           file.path(directory, "grna.odm"), 9876L)
    expect_identical(storage_read_response(reopened$response_matrix, 1L),
                     storage_read_response(pair$response_matrix, 1L))
    expect_error(storage_reopen_backed_pair(file.path(directory, "response.odm"),
                                            file.path(directory, "grna.odm"), 9877L),
                 "distinct IDs", fixed = TRUE)
})
