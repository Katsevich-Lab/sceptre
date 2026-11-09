# data.table::fread() reads a compressed file only through the R.utils
# package, which is not a dependency of sceptre (a user who installs sceptre
# need not have it). The importers must read compressed Cell Ranger files
# without it: here fread() refuses any compressed file path.

refuse_compressed_paths <- function(env = parent.frame()) {
    fread <- data.table::fread
    local_mocked_bindings(
        fread = function(..., file = NULL) {
            if (!is.null(file) && grepl("\\.(gz|bz2|xz)$", file)) {
                stop("fread() would need R.utils to read ", file)
            }
            if (is.null(file)) fread(...) else fread(..., file = file)
        },
        .package = "data.table",
        .env = env
    )
}


test_that("import_data_from_cellranger() reads .gz files without R.utils", {
    directories <- list.files(
        system.file("extdata", "highmoi_example", package = "sceptre"),
        full.names = TRUE
    )
    expect_true(any(grepl(
        "\\.gz$",
        list.files(directories[1])
    )))
    data(grna_target_data_frame_highmoi)
    expected <- import_data_from_cellranger(
        directories = directories,
        grna_target_data_frame = grna_target_data_frame_highmoi,
        moi = "high"
    )
    refuse_compressed_paths()
    for (backed in c(FALSE, TRUE)) {
        sceptre_object <- import_data_from_cellranger(
            directories = directories,
            grna_target_data_frame = grna_target_data_frame_highmoi,
            moi = "high",
            backed = backed,
            directory_to_write = if (backed) withr::local_tempdir()
        )
        expect_equal(
            as.matrix(get_response_matrix(sceptre_object)),
            as.matrix(get_response_matrix(expected))
        )
        expect_equal(
            get_cell_covariates(sceptre_object),
            get_cell_covariates(expected)
        )
    }
})


test_that("read_delimited_file() reads compressed and plain files alike", {
    lines <- c("a\tb\tc", "d\te\tf")
    plain <- withr::local_tempfile(fileext = ".tsv")
    writeLines(lines, plain)
    compressed <- withr::local_tempfile(fileext = ".tsv.gz")
    con <- gzfile(compressed, "w")
    writeLines(lines, con)
    close(con)
    refuse_compressed_paths()
    expect_identical(
        read_delimited_file(compressed, header = FALSE),
        read_delimited_file(plain, header = FALSE)
    )
    empty <- withr::local_tempfile(fileext = ".tsv.gz")
    close(gzfile(empty, "w"))
    expect_error(read_delimited_file(empty, header = FALSE), "is empty")
})
