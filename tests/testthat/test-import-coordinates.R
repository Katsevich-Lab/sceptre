# eminem checks the Matrix Market format itself; these tests cover how sceptre
# builds its matrices from the parsed entries.

write_cellranger_coordinates <- function(directory, entries, n_rows = 5L) {
    dir.create(directory)
    writeLines(
        c(
            "gene-z\tZ\tGene Expression",
            "guide-a\tA\tCRISPR Guide Capture",
            "gene-a\tA\tGene Expression",
            "guide-b\tB\tCRISPR Guide Capture",
            "gene-m\tM\tGene Expression"
        ),
        file.path(directory, "features.tsv")
    )
    writeLines(c("cell-z", "cell-a"), file.path(directory, "barcodes.tsv"))
    # Cell Ranger writes gzipped matrices
    connection <- gzfile(file.path(directory, "matrix.mtx.gz"), "w")
    writeLines(
        c(
            "%%MatrixMarket matrix coordinate real general",
            "% fixture",
            paste(n_rows, 2L, length(entries)),
            entries
        ),
        connection
    )
    close(connection)
    directory
}

cellranger_metadata <- data.frame(
    grna_id = c("guide-a", "guide-b"),
    grna_target = c("target-a", "non-targeting")
)

test_that("Cell Ranger import sorts entries and maps interleaved feature types", {
    directory <- write_cellranger_coordinates(
        file.path(withr::local_tempdir(), "interleaved"),
        c(
            "5 2 3",
            "2 1 8",
            "1 2 0",
            "4 2 4",
            "3 1 2",
            "1 1 5",
            "2 2 1.5",
            "4 1 6",
            "3 2 7"
        )
    )
    object <- suppressMessages(import_data_from_cellranger(
        directory,
        "low",
        cellranger_metadata
    ))
    response <- get_response_matrix(object)
    guides <- get_grna_matrix(object)
    expect_true(methods::validObject(response))
    expect_true(methods::validObject(guides))
    expect_identical(rownames(response), c("gene-z", "gene-a", "gene-m"))
    expect_identical(rownames(guides), c("guide-a", "guide-b"))
    expect_identical(colnames(response), c("cell-z", "cell-a"))
    expect_identical(response@p, c(0L, 2L, 4L, 5L))
    expect_identical(response@j, c(0L, 1L, 0L, 1L, 1L))
    expect_identical(response@x, c(5, 0, 2, 7, 3))
    expect_identical(guides@p, c(0L, 2L, 4L))
    expect_identical(guides@j, c(0L, 1L, 0L, 1L))
    expect_identical(guides@x, c(8, 1.5, 6, 4))

    # The same coordinate in a different directory is a different cell.
    batches <- suppressMessages(import_data_from_cellranger(
        rep(directory, 2L),
        "low",
        cellranger_metadata
    ))
    expect_identical(
        colnames(get_response_matrix(batches)),
        c("b1_cell-z", "b1_cell-a", "b2_cell-z", "b2_cell-a")
    )
})

test_that("Parse import sorts entries and reads cells from rows", {
    directory <- withr::local_tempdir()
    gene <- file.path(directory, "gene.mtx")
    guide <- file.path(directory, "guide.mtx")
    gene_ids <- file.path(directory, "genes.csv")
    guide_ids <- file.path(directory, "guides.csv")
    header <- c(
        "%%MatrixMarket matrix coordinate real general",
        "% fixture",
        "2 2 4"
    )
    writeLines(c(header, "2 2 7", "1 1 5", "2 1 0", "1 2 2"), gene)
    writeLines(c(header, "2 1 1.5", "1 2 6", "2 2 4", "1 1 8"), guide)
    writeLines(
        c("feature_id,feature_name,genome", "gene-z,Z,human", "gene-a,A,human"),
        gene_ids
    )
    writeLines(
        c(
            "feature_id,feature_name,genome",
            "a,guide-a,human",
            "b,guide-b,human"
        ),
        guide_ids
    )
    object <- suppressMessages(import_data_from_parse(
        gene,
        guide,
        gene_ids,
        guide_ids,
        "low",
        cellranger_metadata
    ))
    response <- get_response_matrix(object)
    expect_identical(response@p, c(0L, 2L, 4L))
    expect_identical(response@j, c(0L, 1L, 0L, 1L))
    expect_identical(response@x, c(5, 0, 2, 7))
    expect_identical(rownames(response), c("gene-z", "gene-a"))
    expect_identical(get_grna_matrix(object)@x, c(8, 1.5, 6, 4))
})

test_that("import rejects repeated entries and malformed matrix files", {
    import <- function(entries, n_rows = 5L) {
        directory <- write_cellranger_coordinates(
            file.path(withr::local_tempdir(), "bad"),
            entries,
            n_rows
        )
        suppressMessages(import_data_from_cellranger(
            directory,
            "low",
            cellranger_metadata
        ))
    }
    expect_error(
        import(c("1 1 2", "2 1 1", "1 1 0")),
        "the same (row, column) entry appears more than once",
        fixed = TRUE
    )
    # the parser's own errors reach the user with the file name
    expect_error(import(c("1 1 2", "2.0 1 1")), "matrix.mtx.gz: .*line 5")
    expect_error(
        import(c("1 1 2", "3 2 1"), n_rows = 4L),
        "has 4 rows, but its features.tsv file has 5",
        fixed = TRUE
    )

    # a later features.tsv that repeats the first one does not match it
    first <- write_cellranger_coordinates(
        file.path(withr::local_tempdir(), "first"),
        "1 1 2"
    )
    second <- write_cellranger_coordinates(
        file.path(withr::local_tempdir(), "second"),
        "1 1 2"
    )
    features <- readLines(file.path(second, "features.tsv"))
    writeLines(rep(features, 2L), file.path(second, "features.tsv"))
    expect_error(
        suppressMessages(import_data_from_cellranger(
            c(first, second),
            "low",
            cellranger_metadata
        )),
        "The features.tsv files must match across directories.",
        fixed = TRUE
    )
})
