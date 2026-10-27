test_that("file CSC construction sorts entries without dropping zeros", {
    dt <- data.table::data.table(
        cell_idx = c(1L, 0L, 1L, 0L),
        feature_idx = c(1L, 1L, 0L, 0L),
        x = c(11, 1, 10, 0)
    )
    out <- create_csc_matrix(dt, 3L, 2L)
    expect_s4_class(out, "dgCMatrix")
    expect_true(methods::validObject(out))
    expect_identical(out@Dim, c(2L, 3L))
    expect_identical(out@p, c(0L, 2L, 4L, 4L))
    expect_identical(out@i, c(0L, 1L, 0L, 1L))
    expect_identical(out@x, c(0, 1, 10, 11))

    # Previously this ordering left only one of the three entries visible,
    # although the matrix passed Matrix's validity check.
    missing <- data.table::data.table(
        cell_idx = c(1L, 0L, 1L),
        feature_idx = c(0L, 1L, 2L),
        x = c(2, 3, 4)
    )
    out <- create_csc_matrix(missing, 2L, 3L)
    expect_identical(out@p, c(0L, 1L, 3L))
    expect_identical(out@i, c(1L, 0L, 2L))
    expect_identical(out@x, c(3, 2, 4))
    expect_identical(as.matrix(out), matrix(c(0, 3, 0, 2, 0, 4), 3L, 2L))
})

test_that("file CSC construction checks coordinates and keeps empty files", {
    for (field in c("cell_idx", "feature_idx")) {
        for (bad in c(-1L, 2L, NA_integer_)) {
            invalid <- data.table::data.table(
                cell_idx = 0L,
                feature_idx = 0L,
                x = 1
            )
            invalid[[field]] <- bad
            expect_error(
                create_csc_matrix(invalid, 2L, 2L),
                "outside the matrix dimensions",
                fixed = TRUE
            )
        }
    }
    empty <- data.table::data.table(
        cell_idx = integer(),
        feature_idx = integer(),
        x = numeric()
    )
    out <- create_csc_matrix(empty, 3L, 2L)
    expect_identical(out@p, rep(0L, 4L))
    expect_identical(out@i, integer())
    expect_identical(out@x, numeric())
    expect_identical(out@Dim, c(2L, 3L))
})

test_that("Parse import sorts unsorted file coordinates", {
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
    metadata <- data.frame(
        grna_id = c("guide-a", "guide-b"),
        grna_target = c("target-a", "non-targeting")
    )
    object <- suppressMessages(import_data_from_parse(
        gene,
        guide,
        gene_ids,
        guide_ids,
        "low",
        metadata
    ))
    response <- get_response_matrix(object)
    expect_identical(response@p, c(0L, 2L, 4L))
    expect_identical(response@j, c(0L, 1L, 0L, 1L))
    expect_identical(response@x, c(5, 0, 2, 7))
    expect_identical(rownames(response), c("gene-z", "gene-a"))
    expect_identical(get_grna_matrix(object)@x, c(8, 1.5, 6, 4))
})

test_that("the sortedness check orders by cell, then by feature", {
    expect_true(coordinates_sorted(integer(), integer()))
    expect_true(coordinates_sorted(c(0L, 0L, 1L), c(0L, 2L, 1L)))
    expect_false(coordinates_sorted(c(0L, 1L, 0L), c(0L, 0L, 1L)))
    expect_false(coordinates_sorted(c(0L, 0L), c(2L, 1L)))
})

test_that("file CSC construction rejects repeated coordinates", {
    dt <- data.table::data.table(
        cell_idx = c(1L, 0L, 1L),
        feature_idx = c(0L, 1L, 0L),
        x = c(0, 7, 3)
    )
    expect_error(
        create_csc_matrix(dt, 2L, 2L),
        "Repeated (cell, feature) coordinates",
        fixed = TRUE
    )
    expect_false(any_repeated_coordinates(integer(), integer()))
    expect_false(any_repeated_coordinates(c(0L, 0L, 1L), c(0L, 1L, 1L)))
    expect_true(any_repeated_coordinates(c(0L, 1L, 1L), c(2L, 1L, 1L)))
})

test_that("Parse import rejects repeated file coordinates", {
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
    writeLines(c(header, "1 1 5", "2 2 7", "1 1 0", "1 2 2"), gene)
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
    metadata <- data.frame(
        grna_id = c("guide-a", "guide-b"),
        grna_target = c("target-a", "non-targeting")
    )
    expect_error(
        suppressMessages(import_data_from_parse(
            gene,
            guide,
            gene_ids,
            guide_ids,
            "low",
            metadata
        )),
        "Repeated (cell, feature) coordinates",
        fixed = TRUE
    )
})

write_cellranger_coordinates <- function(directory, entries) {
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
    writeLines(
        c(
            "%%MatrixMarket matrix coordinate real general",
            "% fixture",
            paste(5L, 2L, length(entries)),
            entries
        ),
        file.path(directory, "matrix.mtx")
    )
    directory
}

test_that("Cell Ranger import maps interleaved feature types to the right rows", {
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
    metadata <- data.frame(
        grna_id = c("guide-a", "guide-b"),
        grna_target = c("target-a", "non-targeting")
    )
    object <- suppressMessages(import_data_from_cellranger(
        directory,
        "low",
        metadata
    ))
    response <- get_response_matrix(object)
    guides <- get_grna_matrix(object)
    expect_true(methods::validObject(response))
    expect_true(methods::validObject(guides))
    expect_identical(rownames(response), c("gene-z", "gene-a", "gene-m"))
    expect_identical(rownames(guides), c("guide-a", "guide-b"))
    expect_identical(colnames(response), c("cell-z", "cell-a"))
    expect_identical(colnames(guides), c("cell-z", "cell-a"))
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
        metadata
    ))
    expect_identical(
        colnames(get_response_matrix(batches)),
        c("b1_cell-z", "b1_cell-a", "b2_cell-z", "b2_cell-a")
    )
    duplicate <- write_cellranger_coordinates(
        file.path(withr::local_tempdir(), "duplicate"),
        c("1 1 2", "2 1 1", "1 1 0")
    )
    expect_error(
        suppressMessages(import_data_from_cellranger(
            duplicate,
            "low",
            metadata
        )),
        "Repeated (cell, feature) coordinates",
        fixed = TRUE
    )
})
