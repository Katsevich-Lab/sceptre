test_that("matrix metadata parsing preserves plain and gzip results", {
    directory <- withr::local_tempdir()
    lines <- c(
        "%%MatrixMarket matrix coordinate real general",
        "% fixture",
        "2 3 4",
        "1 1 5"
    )
    plain <- file.path(directory, "matrix.mtx")
    zipped <- file.path(directory, "matrix.mtx.gz")
    writeLines(lines, plain)
    connection <- gzfile(zipped, "w")
    writeLines(lines, connection)
    close(connection)
    before <- getAllConnections()
    expected <- list(
        n_features = 2L,
        n_cells = 3L,
        n_nonzero = 4L,
        n_to_skip = 3L
    )
    expect_identical(get_mtx_metadata(plain), expected)
    expect_identical(get_mtx_metadata(zipped), expected)
    expect_identical(
        get_mtx_metadata(plain, c("n_cells", "n_features", "n_nonzero")),
        list(n_cells = 2L, n_features = 3L, n_nonzero = 4L, n_to_skip = 3L)
    )
    expect_identical(getAllConnections(), before)
})

test_that("a matrix file without a dimension line leaves no open connection", {
    file <- withr::local_tempfile()
    writeLines("%%MatrixMarket matrix coordinate real general", file)
    before <- getAllConnections()
    expect_error(get_mtx_metadata(file))
    expect_identical(getAllConnections(), before)
})

test_that("metadata read failure closes only the helper-owned connection", {
    file <- withr::local_tempfile()
    writeLines("%%MatrixMarket matrix coordinate real general", file)
    helper <- get_mtx_metadata
    scope <- new.env(parent = environment(helper))
    captured <- NULL
    scope$readLines <- function(con, ...) {
        captured <<- con
        stop("Controlled metadata read failure.")
    }
    environment(helper) <- scope
    before <- getAllConnections()
    expect_error(
        helper(file),
        "Controlled metadata read failure.",
        fixed = TRUE
    )
    expect_false(as.integer(captured) %in% getAllConnections())
    expect_identical(getAllConnections(), before)
})

test_that("pipeline vector writes preserve lines and close a failed write", {
    path <- withr::local_tempfile()
    writer <- write_vector
    before <- getAllConnections()
    expect_null(writer(c("z", "a", "z"), path))
    expect_identical(readLines(path), c("z", "a", "z"))
    expect_identical(getAllConnections(), before)

    scope <- new.env(parent = environment(writer))
    captured <- NULL
    scope$writeLines <- function(text, con, ...) {
        captured <<- con
        stop("Controlled vector write failure.")
    }
    environment(writer) <- scope
    expect_error(
        writer(1:3, path),
        "Controlled vector write failure.",
        fixed = TRUE
    )
    expect_false(as.integer(captured) %in% getAllConnections())
    expect_identical(getAllConnections(), before)
})
