test_that("triplet conversion leaves the caller's matrix unchanged", {
    coordinates <- list(i = c(2L, 0L, 1L, 0L, 2L), j = c(0L, 2L, 0L, 0L, 2L))
    values <- c(1.25, 0, 4.5, 2, 0.25)
    names <- list(c("r1", "r2", "r3"), c("c1", "c2", "c3"))
    numeric <- methods::new(
        "dgTMatrix",
        i = coordinates$i,
        j = coordinates$j,
        x = values,
        Dim = c(3L, 3L),
        Dimnames = names
    )
    logical <- methods::new(
        "lgTMatrix",
        i = coordinates$i,
        j = coordinates$j,
        x = c(TRUE, FALSE, TRUE, FALSE, TRUE),
        Dim = c(3L, 3L),
        Dimnames = names
    )
    for (input in list(numeric, logical)) {
        before <- serialize(input, NULL)
        aliases <- list(i = input@i, j = input@j, x = input@x)
        aliases_before <- serialize(aliases, NULL)
        for (row_access in c(TRUE, FALSE)) {
            actual <- set_matrix_accessibility(input, row_access)
            order <- if (row_access) {
                order(coordinates$i, coordinates$j)
            } else {
                order(coordinates$j, coordinates$i)
            }
            leading <- if (row_access) {
                coordinates$i[order]
            } else {
                coordinates$j[order]
            }
            counts <- tabulate(leading + 1L, nbins = 3L)
            expect_identical(actual@p, as.integer(c(0, cumsum(counts))))
            expect_identical(
                if (row_access) actual@j else actual@i,
                if (row_access) coordinates$j[order] else coordinates$i[order]
            )
            expect_identical(
                actual@x,
                if (methods::is(input, "lgTMatrix")) {
                    rep(1, 5L)
                } else {
                    values[order]
                }
            )
            expect_identical(
                class(actual)[[1L]],
                if (row_access) "dgRMatrix" else "dgCMatrix"
            )
            expect_identical(dimnames(actual), names)
            expect_identical(serialize(input, NULL), before)
            expect_identical(serialize(aliases, NULL), aliases_before)
        }
    }
})

test_that("import_data() leaves a caller's triplet matrix unchanged", {
    response <- methods::new(
        "dgTMatrix",
        i = c(1L, 0L, 1L, 0L),
        j = c(1L, 1L, 0L, 0L),
        x = c(4, 3, 2, 1),
        Dim = c(2L, 2L),
        Dimnames = list(c("r1", "r2"), c("c1", "c2"))
    )
    guides <- methods::new(
        "dgTMatrix",
        i = c(1L, 0L),
        j = c(1L, 0L),
        x = c(5, 6),
        Dim = c(2L, 2L),
        Dimnames = list(c("g1", "g2"), c("c1", "c2"))
    )
    targets <- data.frame(
        grna_id = c("g1", "g2"),
        grna_target = c("non-targeting", "t")
    )
    before <- serialize(list(response, guides), NULL)
    object <- import_data(response, guides, targets, "low")
    expect_identical(serialize(list(response, guides), NULL), before)
    expect_identical(
        as.matrix(get_response_matrix(object)),
        as.matrix(response)
    )
})
