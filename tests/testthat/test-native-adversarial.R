test_that("sparse native QC agrees with structural membership on adversarial views", {
    withr::with_seed(23702L, {
        for (iteration in seq_len(12L)) {
            dense <- matrix(sample(c(0, 0, 0, 1, 3), 4L * 9L, replace = TRUE), 4L)
            store <- methods::as(Matrix::Matrix(dense, sparse = TRUE), "RsparseMatrix")
            # Native QC deliberately counts stored entries, including explicit
            # zeros. Use structural membership rather than a dense-value oracle.
            if (length(store@x)) store@x[1L] <- 0
            structure <- matrix(FALSE, 4L, 9L)
            for (row in seq_len(4L)) {
                offsets <- seq.int(store@p[row] + 1L, store@p[row + 1L])
                if (store@p[row] < store@p[row + 1L])
                    structure[row, store@j[offsets] + 1L] <- TRUE
            }
            cells <- sample(seq_len(9L), 11L, replace = TRUE)
            groups <- list(a = c(9L, 3L, 6L), b = c(11L, 2L), empty = integer())
            nt <- list(n1 = c(1L, 3L, 5L), n2 = c(3L, 7L, 11L))
            pairs <- data.frame(response_idx = c(4L, 1L, 4L, 2L, 3L),
                                group_idx = c(2L, 1L, 2L, 3L, 1L))
            snapshot <- serialize(list(store, cells, groups, nt, pairs), NULL)
            expected_nt <- vapply(seq_len(4L), function(row)
                vapply(nt, function(group) sum(structure[row, cells[group]]), integer(1)),
                integer(2L))
            expected_trt <- vapply(seq_len(nrow(pairs)), function(k)
                sum(structure[pairs$response_idx[k], cells[groups[[pairs$group_idx[k]]]]]), integer(1))
            for (control in list("complement", unique(unlist(nt)))) {
                result <- sceptre:::storage_pair_counts(store, pairs, groups, cells, control, nt)
                totals <- if (is.character(control)) as.integer(rowSums(structure[, cells, drop = FALSE])) else
                    as.integer(colSums(expected_nt))
                controls <- totals[pairs$response_idx] - if (is.character(control)) expected_trt else 0L
                expect_identical(result$n_nonzero_mat, unname(expected_nt))
                expect_identical(result$n_nonzero_tot, totals)
                expect_identical(result$n_nonzero_trt, expected_trt)
                expect_identical(result$n_nonzero_cntrl, controls)
                expect_identical(serialize(list(store, cells, groups, nt, pairs), NULL), snapshot)
            }
        }
    })
})

test_that("native mitochondrial coordinates reject missing and out-of-range rows", {
    for (row in c(0L, -1L, 3L, NA_integer_)) {
        expect_error(sceptre:::compute_cell_covariates_cpp(0L, c(0L, 1L), 1,
            2L, 1L, row, TRUE, TRUE), "one-based feature range")
    }
})
