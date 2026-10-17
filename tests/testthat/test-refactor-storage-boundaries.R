test_that("count boundaries materialize borrowed compact and coercible coordinates", {
    # R can store sequences lazily (ALTREP). Check these position vectors too,
    # and keep the C++ helper's existing conversion of numeric groups to integers.
    j <- c(0L, 2L, 1L, 3L)
    p <- c(0L, 2L, 4L)
    cells <- 1:4
    groups <- list(all = 1:4, even = c(4, 2), empty = integer())
    nt <- list(1:2, 2:4)
    before <- serialize(list(j, p, cells, groups, nt), NULL)
    expected <- matrix(c(2L, 0L, 0L, 2L, 2L, 0L), 3L, 2L)
    expect_identical(sparse_group_counts(j, p, 4L, groups, cells), expected)
    actual <- sparse_pair_counts(j, p, 4L, groups, nt, 1:2, 2:1, FALSE, cells)
    expect_identical(actual, list(n_nonzero_mat = matrix(c(1L, 1L, 1L, 2L), 2L, 2L),
        n_nonzero_tot = c(2L, 3L), n_nonzero_trt = c(0L, 2L),
        n_nonzero_cntrl = c(2L, 3L)))
    expect_identical(serialize(list(j, p, cells, groups, nt), NULL), before)
    expect_error(sparse_group_counts(j, p, 4L, list(c(1L, NA_integer_)), cells),
                 "one-based coordinate view", fixed = TRUE)
})

test_that("mitochondrial boundary materializes rows before its owned compute scope", {
    rows <- 1:4
    mt <- 1:4
    before <- serialize(list(rows, mt), NULL)
    actual <- compute_cell_covariates_cpp(rows - 1L, c(0L, 4L), as.numeric(rows),
        4L, 1L, mt, TRUE, TRUE)
    expect_identical(actual, list(n_nonzero = 4L, n_umi = 10, p_mito = 1,
        max_feature = 3L, frac_umis_max_feature = 0.4))
    expect_identical(serialize(list(rows, mt), NULL), before)
    disabled <- compute_cell_covariates_cpp(integer(), c(0L, 0L), numeric(),
        4L, 1L, NA_integer_, FALSE, FALSE)
    expect_identical(disabled, list(n_nonzero = 0L, n_umi = 0, p_mito = 0,
        max_feature = 0L, frac_umis_max_feature = 0))
})
