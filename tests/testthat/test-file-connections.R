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
