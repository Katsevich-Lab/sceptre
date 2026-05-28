# First-signal tests for the public importers.
# Loads bundled fixtures shipped in inst/extdata/ and asserts basic shape of
# the returned `sceptre_object`. Exact values are *not* asserted -- only that
# the object is well-formed and that the response/grna matrices have non-zero,
# consistent dimensions. Intended to catch gross breakage in the import path
# without becoming brittle to fixture content.

test_that("import_data_from_cellranger() returns a well-shaped sceptre_object", {
  directories <- c(
    system.file("extdata", "highmoi_example", "gem_group_1", package = "sceptre"),
    system.file("extdata", "highmoi_example", "gem_group_2", package = "sceptre")
  )
  # sanity: fixtures must be present in the installed package
  expect_true(all(nzchar(directories)))
  expect_true(all(dir.exists(directories)))

  data("grna_target_data_frame_highmoi")

  sceptre_object <- suppressMessages(
    import_data_from_cellranger(
      directories = directories,
      moi = "high",
      grna_target_data_frame = grna_target_data_frame_highmoi
    )
  )

  # type
  expect_s4_class(sceptre_object, "sceptre_object")

  # matrices: both modalities present, non-empty, and column counts agree
  response_matrix <- get_response_matrix(sceptre_object)
  grna_matrix <- get_grna_matrix(sceptre_object)
  expect_true(nrow(response_matrix) > 0L)
  expect_true(ncol(response_matrix) > 0L)
  expect_true(nrow(grna_matrix) > 0L)
  expect_equal(ncol(response_matrix), ncol(grna_matrix))

  # covariate data frame: one row per cell
  covariate_df <- get_cell_covariates(sceptre_object)
  expect_s3_class(covariate_df, "data.frame")
  expect_equal(nrow(covariate_df), ncol(response_matrix))

  # batch column is populated when there are >= 2 directories
  expect_true("batch" %in% colnames(covariate_df))
  expect_equal(length(unique(covariate_df$batch)), length(directories))

  # grna_target_data_frame round-trips into the object
  expect_s3_class(sceptre_object@grna_target_data_frame, "data.frame")
  expect_true(nrow(sceptre_object@grna_target_data_frame) > 0L)

  # MOI flag wired through
  expect_false(sceptre_object@low_moi)
})


test_that("import_data_from_parse() returns a well-shaped sceptre_object", {
  directory <- system.file("extdata", "parse_example", package = "sceptre")
  expect_true(nzchar(directory))

  gene_mat_fp <- file.path(directory, "gene_mat.mtx")
  grna_mat_fp <- file.path(directory, "grna_mat.mtx")
  all_genes_fp <- file.path(directory, "all_genes.csv")
  all_grnas_fp <- file.path(directory, "all_grnas.csv")
  expect_true(all(file.exists(c(gene_mat_fp, grna_mat_fp, all_genes_fp, all_grnas_fp))))

  # Matches the @examples block of import_data_from_parse().
  grna_target_data_frame <- data.frame(
    grna_id = c("guide_A", "guide_B", "guide_C"),
    grna_target = c("target-A", "target-B", "non-targeting")
  )

  sceptre_object <- suppressMessages(
    import_data_from_parse(
      gene_mat_fp = gene_mat_fp,
      grna_mat_fp = grna_mat_fp,
      all_genes_fp = all_genes_fp,
      all_grnas_fp = all_grnas_fp,
      moi = "low",
      grna_target_data_frame = grna_target_data_frame
    )
  )

  # type
  expect_s4_class(sceptre_object, "sceptre_object")

  # matrices: both modalities present, non-empty, and column counts agree
  response_matrix <- get_response_matrix(sceptre_object)
  grna_matrix <- get_grna_matrix(sceptre_object)
  expect_true(nrow(response_matrix) > 0L)
  expect_true(ncol(response_matrix) > 0L)
  expect_true(nrow(grna_matrix) > 0L)
  expect_equal(ncol(response_matrix), ncol(grna_matrix))

  # grna matrix row count matches the supplied target data frame
  expect_equal(nrow(grna_matrix), nrow(grna_target_data_frame))

  # covariate data frame: one row per cell
  covariate_df <- get_cell_covariates(sceptre_object)
  expect_s3_class(covariate_df, "data.frame")
  expect_equal(nrow(covariate_df), ncol(response_matrix))

  # MOI flag wired through
  expect_true(sceptre_object@low_moi)
})
