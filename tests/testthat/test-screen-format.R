# screen.h5mu carries a format stamp in /uns/sceptre_format, a JSON object
# whose `format` names the layout. sceptre 0.99.2 writes format 1, the first
# public layout (registry: contracts/file-formats.md in sceptre-hub). A reader
# accepts every format up to its own and refuses a later one.

test_that("screen.h5mu is stamped as format 1", {
    data(lowmoi_example_data, package = "sceptre", envir = environment())
    sceptre_object <- import_data(
        response_matrix = lowmoi_example_data$response_matrix,
        grna_matrix = lowmoi_example_data$grna_matrix,
        grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
        extra_covariates = lowmoi_example_data$extra_covariates,
        moi = "low"
    )
    directory <- withr::local_tempdir()
    write_sceptre_object(sceptre_object, directory)
    stamp <- rhdf5::h5read(file.path(directory, "screen.h5mu"), "/uns/sceptre_format")
    expect_identical(as.character(stamp), "{\"format\":1}")
    expect_identical(sceptre_file_format, 1L)
})


test_that("the stamp is read as JSON and later formats are refused", {
    expect_silent(check_screen_format("{\"format\":1}", "f"))
    expect_silent(check_screen_format("{ \"format\" : 1 }", "f"))
    expect_silent(check_screen_format("{\"format\":1,\"readers\":{}}", "f"))
    expect_error(check_screen_format("{\"format\":2}", "f"), "newer version")
    expect_error(check_screen_format("{\"format\":12}", "f"), "file format 12")
    expect_error(check_screen_format("{\"version\":1}", "f"), "not a screen.h5mu")
    expect_error(check_screen_format(NULL, "f"), "not a screen.h5mu")
    expect_error(check_screen_format(c("{\"format\":1}", "x"), "f"), "not a screen.h5mu")
})
