get_cellranger_test_inputs <- function() {
    data("grna_target_data_frame_highmoi")
    directories <- paste0(
        system.file("extdata", package = "sceptre"),
        "/highmoi_example/gem_group_",
        c(1, 2)
    )
    list(
        directories = directories,
        grna_target_data_frame = grna_target_data_frame_highmoi
    )
}


copy_cellranger_input_files <- function(source_directory, output_directory) {
    dir.create(output_directory)
    input_files <- list.files(
        source_directory,
        pattern = "^(features\\.tsv|matrix\\.mtx)(\\.gz)?$",
        full.names = TRUE
    )
    file.copy(input_files, output_directory)
}


test_that("import_data_from_cellranger retains cell barcodes", {
    inputs <- get_cellranger_test_inputs()
    expected_barcodes <- unlist(
        lapply(inputs$directories, function(directory) {
            data.table::fread(
                file.path(directory, "barcodes.tsv.gz"),
                header = FALSE,
                colClasses = "character"
            )[[1L]]
        }),
        use.names = FALSE
    )

    sceptre_object <- suppressMessages(import_data_from_cellranger(
        directories = inputs$directories,
        moi = "high",
        grna_target_data_frame = inputs$grna_target_data_frame
    ))

    expect_identical(
        colnames(get_response_matrix(sceptre_object)),
        expected_barcodes
    )
    expect_identical(
        colnames(get_grna_matrix(sceptre_object)),
        expected_barcodes
    )
    expect_identical(
        rownames(get_cell_covariates(sceptre_object)),
        expected_barcodes
    )
})


test_that("Cell Ranger barcode files remain optional", {
    inputs <- get_cellranger_test_inputs()
    directory <- file.path(tempdir(), "cellranger_without_barcodes")
    unlink(directory, recursive = TRUE)
    copy_cellranger_input_files(inputs$directories[1L], directory)

    sceptre_object <- suppressMessages(import_data_from_cellranger(
        directories = directory,
        moi = "high",
        grna_target_data_frame = inputs$grna_target_data_frame
    ))

    expect_null(colnames(get_response_matrix(sceptre_object)))
    expect_null(colnames(get_grna_matrix(sceptre_object)))
})


test_that("Cell Ranger barcode files cannot be supplied only partially", {
    inputs <- get_cellranger_test_inputs()
    directory <- file.path(tempdir(), "cellranger_partial_barcodes")
    unlink(directory, recursive = TRUE)
    copy_cellranger_input_files(inputs$directories[1L], directory)

    expect_error(
        suppressMessages(import_data_from_cellranger(
            directories = c(directory, inputs$directories[2L]),
            moi = "high",
            grna_target_data_frame = inputs$grna_target_data_frame
        )),
        paste0(
            "Barcode files must either be present in every input directory ",
            "or absent from every input directory"
        ),
        fixed = TRUE
    )
})


test_that("duplicate barcodes across Cell Ranger directories are disambiguated", {
    inputs <- get_cellranger_test_inputs()
    directory <- inputs$directories[1L]
    barcodes <- data.table::fread(
        file.path(directory, "barcodes.tsv.gz"),
        header = FALSE,
        colClasses = "character"
    )[[1L]]
    expected_barcodes <- c(
        paste0("b1_", barcodes),
        paste0("b2_", barcodes)
    )

    sceptre_object <- suppressMessages(import_data_from_cellranger(
        directories = rep(directory, 2L),
        moi = "high",
        grna_target_data_frame = inputs$grna_target_data_frame
    ))

    expect_identical(
        colnames(get_response_matrix(sceptre_object)),
        expected_barcodes
    )
    expect_identical(
        colnames(get_grna_matrix(sceptre_object)),
        expected_barcodes
    )
    expect_identical(
        rownames(get_cell_covariates(sceptre_object)),
        expected_barcodes
    )
})


test_that("Cell Ranger barcodes agree with named extra covariates", {
    inputs <- get_cellranger_test_inputs()
    directory <- inputs$directories[1L]
    barcodes <- data.table::fread(
        file.path(directory, "barcodes.tsv.gz"),
        header = FALSE,
        colClasses = "character"
    )[[1L]]
    extra_covariates <- data.frame(covariate = seq_along(barcodes))
    rownames(extra_covariates) <- barcodes

    sceptre_object <- suppressMessages(import_data_from_cellranger(
        directories = directory,
        moi = "high",
        grna_target_data_frame = inputs$grna_target_data_frame,
        extra_covariates = extra_covariates
    ))
    expect_identical(
        rownames(get_cell_covariates(sceptre_object)),
        barcodes
    )

    rownames(extra_covariates) <- paste0("mismatch_", barcodes)
    expect_error(
        suppressMessages(import_data_from_cellranger(
            directories = directory,
            moi = "high",
            grna_target_data_frame = inputs$grna_target_data_frame,
            extra_covariates = extra_covariates
        )),
        paste0(
            "provided cell barcodes in the `response_matrix` and ",
            "`extra_covariates`"
        ),
        fixed = TRUE
    )
})


test_that("Cell Ranger barcode counts must match matrix cell counts", {
    inputs <- get_cellranger_test_inputs()
    directory <- file.path(tempdir(), "cellranger_bad_barcodes")
    unlink(directory, recursive = TRUE)
    copy_cellranger_input_files(inputs$directories[1L], directory)
    writeLines("cell_1", file.path(directory, "barcodes.tsv"))

    expect_error(
        suppressMessages(import_data_from_cellranger(
            directories = directory,
            moi = "high",
            grna_target_data_frame = inputs$grna_target_data_frame
        )),
        "contains 1 barcode, but its matrix contains 242 cells",
        fixed = TRUE
    )
})
