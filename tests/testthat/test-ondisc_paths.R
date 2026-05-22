test_that("ondisc output directories are normalized before use", {
  directory_to_write <- file.path(tempdir(), "sceptre ondisc path")
  prepared_directory <- prepare_directory_to_write(directory_to_write)

  expect_true(dir.exists(prepared_directory))
  expect_false(grepl("\\\\", prepared_directory))
  expect_false(grepl("~", prepared_directory, fixed = TRUE))
})

test_that("ondisc output paths reject tildes after normalization", {
  directory_to_write <- file.path(tempdir(), "sceptre~ondisc")
  dir.create(directory_to_write, showWarnings = FALSE)

  expect_error(
    prepare_directory_to_write(directory_to_write),
    "contains a '~' character inside a directory or file name",
    fixed = TRUE
  )
})

test_that("leading ~/ paths are accepted (expanded before the tilde check)", {
  directory_to_write <- "~/sceptre_leading_tilde_test"
  on.exit(unlink(path.expand(directory_to_write), recursive = TRUE), add = TRUE)

  prepared <- prepare_directory_to_write(directory_to_write)

  expect_true(dir.exists(prepared))
  expect_false(grepl("~", prepared, fixed = TRUE))
})

test_that("import_data can create ondisc-backed objects in a normalized directory", {
  testthat::skip_if_not_installed("ondisc")
  data("lowmoi_example_data")
  directory_to_write <- file.path(tempdir(), "sceptre ondisc import")

  sceptre_object <- import_data(
    response_matrix = lowmoi_example_data$response_matrix,
    grna_matrix = lowmoi_example_data$grna_matrix,
    grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
    extra_covariates = lowmoi_example_data$extra_covariates,
    moi = "low",
    use_ondisc = TRUE,
    directory_to_write = directory_to_write
  )

  prepared_directory <- prepare_directory_to_write(directory_to_write)
  expect_true(file.exists(file.path(prepared_directory, "response.odm")))
  expect_true(file.exists(file.path(prepared_directory, "grna.odm")))
  expect_true(methods::is(get_response_matrix(sceptre_object), "odm"))
  expect_true(methods::is(get_grna_matrix(sceptre_object), "odm"))
})
