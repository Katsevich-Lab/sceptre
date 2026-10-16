test_that("canonical counts preserve repeated source and pair views without mutation", {
    store <- storage_test_matrix()
    # Independently spell out stored source columns, including stored zero values.
    stored <- list(c(1L, 4L, 8L, 11L), c(2L, 5L, 8L, 12L))
    cells <- c(11L, 4L, 8L, 2L, 11L)
    groups <- list(A = c(2L, 4L), B = c(5L, 3L, 1L), empty = integer())
    nt <- list(nt1 = c(1L, 2L), nt2 = c(2L, 3L))
    pairs <- data.frame(response_idx = c(2L, 1L, 2L, 1L, 1L),
                        group_idx = c(2L, 1L, 1L, 2L, 3L))
    before <- serialize(list(store, cells, groups, nt, pairs), NULL)
    group_matrix <- vapply(stored, function(present) {
        vapply(nt, function(group) sum(cells[group] %in% present), integer(1))
    }, integer(length(nt)))
    treated <- vapply(seq_len(nrow(pairs)), function(row) {
        sum(cells[groups[[pairs$group_idx[[row]]]]] %in% stored[[pairs$response_idx[[row]]]])
    }, integer(1))
    for (controls in list("complement", unique(unlist(nt)))) {
        result <- storage_pair_counts(store, pairs, groups, cells, controls, nt)
        totals <- if (is.character(controls)) vapply(stored, function(present) {
            sum(cells %in% present)
        }, integer(1)) else as.integer(colSums(group_matrix))
        expect_identical(result$n_nonzero_trt, treated)
        expect_identical(result$n_nonzero_mat, unname(group_matrix))
        expect_identical(result$n_nonzero_tot, unname(totals))
        expected_controls <- totals[pairs$response_idx]
        if (is.character(controls)) expected_controls <- expected_controls - treated
        expect_identical(result$n_nonzero_cntrl, unname(expected_controls))
        expect_identical(serialize(list(store, cells, groups, nt, pairs), NULL), before)
    }
    # NT totals deliberately sum overlapping groups: analysis position 2 occurs
    # twice, unlike a union count. The native kernel must not normalize this.
    summed <- storage_pair_counts(store, pairs, groups, cells, 1:3, nt)
    expect_identical(summed$n_nonzero_tot[[1L]], 4L)
    expect_identical(sum(cells[unique(unlist(nt))] %in% stored[[1L]]), 3L)
    expected <- vapply(stored, function(present) {
        vapply(groups, function(group) sum(cells[group] %in% present), integer(1))
    }, integer(length(groups)))
    for (requested in list(c(2L, 1L, 2L), 1L, integer())) {
        expect_equal(storage_group_counts(store, requested, groups, cells),
                     unname(expected[, requested, drop = FALSE]))
    }
    expect_identical(serialize(list(store, cells, groups, nt, pairs), NULL), before)
})

test_that("empty sparse count views preserve result dimensions", {
    store <- storage_test_matrix()
    empty_pairs <- data.frame(response_idx = integer(), group_idx = integer())
    out <- storage_pair_counts(store, empty_pairs, list(), integer(), "complement", list())
    expect_identical(dim(out$n_nonzero_mat), c(0L, 2L))
    expect_identical(out$n_nonzero_tot, c(0L, 0L))
    expect_identical(out$n_nonzero_trt, integer())
    expect_identical(out$n_nonzero_cntrl, integer())
    expect_identical(dim(storage_group_counts(store, c(2L, 1L), list(), integer())), c(0L, 2L))
    expect_identical(storage_group_counts(store, 1L, list(integer()), integer()), matrix(0L, 1L, 1L))
})

test_that("mitochondrial coordinates are owned and value summaries remain distinct", {
    store <- set_matrix_accessibility(storage_test_matrix(), FALSE)
    mt <- c(2L, 1L, 2L)
    before <- serialize(list(store, mt), NULL)
    loaded <- as.matrix(store)
    actual <- compute_cell_covariates_cpp(store@i, store@p, store@x,
        nrow(store), ncol(store), mt, TRUE, TRUE)
    expect_equal(actual$n_umi, colSums(loaded))
    expect_identical(actual$n_nonzero, as.integer(colSums(loaded > 0.5)))
    expect_equal(actual$p_mito, colSums(loaded[mt, , drop = FALSE]) / colSums(loaded))
    expect_identical(actual$max_feature, apply(loaded, 2L, which.max) - 1L)
    expect_equal(actual$frac_umis_max_feature, apply(loaded, 2L, max) / colSums(loaded))
    expect_identical(storage_cell_covariates(store, mt, TRUE, TRUE), actual)
    expect_identical(compute_cell_covariates_cpp(store@i, store@p, store@x,
        nrow(store), ncol(store), mt, TRUE, TRUE), actual)
    expect_identical(serialize(list(store, mt), NULL), before)
})

test_that("triplet conversion owns sorting and preserves stored entry semantics", {
    coordinates <- list(i = c(2L, 0L, 1L, 0L, 2L), j = c(0L, 2L, 0L, 0L, 2L))
    values <- c(1.25, 0, 4.5, 2, 0.25)
    names <- list(c("r1", "r2", "r3"), c("c1", "c2", "c3"))
    numeric <- methods::new("dgTMatrix", i = coordinates$i, j = coordinates$j,
        x = values, Dim = c(3L, 3L), Dimnames = names)
    logical <- methods::new("lgTMatrix", i = coordinates$i, j = coordinates$j,
        x = c(TRUE, FALSE, TRUE, FALSE, TRUE), Dim = c(3L, 3L), Dimnames = names)
    repeated <- methods::new("dgTMatrix", i = coordinates$i,
        j = c(0L, 2L, 0L, 0L, 0L), x = values, Dim = c(3L, 3L), Dimnames = names)
    for (input in list(numeric, logical, repeated)) {
        coordinates <- list(i = input@i, j = input@j)
        before <- serialize(input, NULL)
        aliases <- list(i = input@i, j = input@j, x = input@x)
        aliases_before <- serialize(aliases, NULL)
        for (row_access in c(TRUE, FALSE)) {
            actual <- set_matrix_accessibility(input, row_access)
            order <- if (row_access) order(coordinates$i, coordinates$j) else
                order(coordinates$j, coordinates$i)
            leading <- if (row_access) coordinates$i[order] else coordinates$j[order]
            counts <- tabulate(leading + 1L, nbins = 3L)
            expect_identical(actual@p, as.integer(c(0, cumsum(counts))))
            expect_identical(if (row_access) actual@j else actual@i,
                             if (row_access) coordinates$j[order] else coordinates$i[order])
            expect_identical(actual@x, if (methods::is(input, "lgTMatrix")) rep(1, 5L) else values[order])
            expect_identical(class(actual)[[1L]], if (row_access) "dgRMatrix" else "dgCMatrix")
            expect_identical(dimnames(actual), names)
            expect_identical(serialize(input, NULL), before)
            expect_identical(serialize(aliases, NULL), aliases_before)
        }
    }
    expect_identical(sparse_pointer_to_indices(c(0L, 0L, 2L, 2L, 3L)), c(1L, 1L, 3L))
})

test_that("sorting preserves borrowed value bits and every source representation", {
    # A nondefault quiet-NaN payload and negative zero must survive duplication;
    # arithmetic used as a copying shortcut could change these bits.
    payload <- readBin(as.raw(c(0x17, 0, 0, 0, 0, 0, 0xf8, 0x7f)),
                       what = "double", n = 1L, size = 8L, endian = "little")
    values <- c(-0, payload, 1.25, -Inf, Inf)
    triplet <- methods::new("dgTMatrix", i = c(2L, 0L, 1L, 0L, 2L),
        j = c(0L, 2L, 0L, 0L, 2L), x = values, Dim = c(3L, 3L))
    before <- serialize(triplet, NULL)
    for (row_access in c(TRUE, FALSE)) {
        order <- if (row_access) order(triplet@i, triplet@j) else order(triplet@j, triplet@i)
        actual <- set_matrix_accessibility(triplet, row_access)
        expect_identical(writeBin(actual@x, raw(), size = 8L),
                         writeBin(values[order], raw(), size = 8L))
        expect_identical(serialize(triplet, NULL), before)
    }

    dense <- matrix(c(3, 1, 4, 2, 6, 5), 2L, 3L,
                    dimnames = list(c("r1", "r2"), c("c1", "c2", "c3")))
    numeric_inputs <- list(dense, methods::as(dense, "RsparseMatrix"),
                          methods::as(dense, "CsparseMatrix"))
    for (input in numeric_inputs) {
        before <- serialize(input, NULL)
        aliases <- if (methods::is(input, "sparseMatrix")) {
            list(index = if (methods::is(input, "RsparseMatrix")) input@j else input@i,
                 pointer = input@p, value = input@x)
        } else list(value = input)
        aliases_before <- serialize(aliases, NULL)
        for (row_access in c(TRUE, FALSE)) {
            actual <- set_matrix_accessibility(input, row_access)
            expect_identical(as.matrix(actual), dense)
            expect_identical(serialize(input, NULL), before)
            expect_identical(serialize(aliases, NULL), aliases_before)
        }
    }

    names <- list(c("r1", "r2", "r3"), c("c1", "c2", "c3"))
    logical_csc <- methods::new("lgCMatrix", i = c(1L, 0L, 2L),
        p = 0:3, x = c(FALSE, TRUE, FALSE), Dim = c(3L, 3L), Dimnames = names)
    logical_csr <- methods::as(logical_csc, "RsparseMatrix")
    expected <- matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 1), 3L, 3L, dimnames = names)
    for (input in list(logical_csc, logical_csr)) {
        before <- serialize(input, NULL)
        for (row_access in c(TRUE, FALSE)) {
            actual <- set_matrix_accessibility(input, row_access)
            expect_identical(as.matrix(actual), expected)
            expect_identical(serialize(input, NULL), before)
        }
    }
})
