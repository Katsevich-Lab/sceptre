test_that("dense conversion preserves symmetric and near-symmetric entries", {
    symmetric <- matrix(
        c(1, 2, 2, 3),
        2L,
        2L,
        dimnames = list(c("c1", "c2"), c("c1", "c2"))
    )
    near <- symmetric
    near[2L, 1L] <- 2 + 4 * .Machine$double.eps
    for (input in list(symmetric, near)) {
        before <- serialize(input, NULL)
        for (row_access in c(TRUE, FALSE)) {
            actual <- set_matrix_accessibility(input, row_access)
            expect_identical(as.matrix(actual), input)
            expect_identical(
                actual@x,
                if (row_access) as.vector(t(input)) else as.vector(input)
            )
            expect_identical(
                class(actual)[[1L]],
                if (row_access) "dgRMatrix" else "dgCMatrix"
            )
            expect_identical(serialize(input, NULL), before)
        }
    }
    # non-square dense input is converted directly, as before
    wide <- matrix(c(0, 2, 3, 0, 0, 5), 2L, 3L)
    actual <- set_matrix_accessibility(wide, TRUE)
    expect_identical(as.matrix(actual), wide)
    expect_identical(class(actual)[[1L]], "dgRMatrix")
})

test_that("dense import retains both triangles and correct cell totals", {
    response <- matrix(
        c(1, 2, 2, 3),
        2L,
        2L,
        dimnames = list(c("c1", "c2"), c("c1", "c2"))
    )
    guides <- matrix(
        c(1, 0, 0, 1),
        2L,
        2L,
        dimnames = list(c("g1", "g2"), c("c1", "c2"))
    )
    targets <- data.frame(
        grna_id = c("g1", "g2"),
        grna_target = c("non-targeting", "target")
    )
    for (lower in c(2, 2 + 4 * .Machine$double.eps)) {
        response[2L, 1L] <- lower
        before <- serialize(list(response, guides, targets), NULL)
        object <- import_data(
            response,
            guides,
            targets,
            "low",
            response_names = c("MT-first", "second")
        )
        covariates <- get_cell_covariates(object)
        expect_identical(as.matrix(get_response_matrix(object)), response)
        expect_identical(covariates$response_n_umis, unname(colSums(response)))
        expect_identical(covariates$response_n_nonzero, c(2L, 2L))
        expect_equal(
            covariates$response_p_mito,
            unname(response[1L, ] / colSums(response))
        )
        expect_identical(
            serialize(list(response, guides, targets), NULL),
            before
        )
    }
})
