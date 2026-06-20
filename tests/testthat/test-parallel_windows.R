test_that("parallel mode is rejected on Windows", {
    expect_error(
        check_parallel_supported(TRUE, os_type = "windows"),
        "`parallel = TRUE` is not supported on Windows",
        fixed = TRUE
    )
    expect_no_error(check_parallel_supported(FALSE, os_type = "windows"))
    expect_no_error(check_parallel_supported(TRUE, os_type = "unix"))
})

test_that("parallel argument must be a single logical value", {
    expect_error(
        check_parallel_supported(NA),
        "`parallel` should be a single logical value.",
        fixed = TRUE
    )
    expect_error(
        check_parallel_supported(c(TRUE, FALSE)),
        "`parallel` should be a single logical value.",
        fixed = TRUE
    )
    expect_error(
        check_parallel_supported("TRUE"),
        "`parallel` should be a single logical value.",
        fixed = TRUE
    )
})
