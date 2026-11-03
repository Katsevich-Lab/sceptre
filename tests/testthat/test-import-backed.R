# A small Cell Ranger directory: features.tsv with 3 genes and 2 gRNAs (gene
# and gRNA rows interleaved), and matrix.mtx holding `entries` (feature, cell,
# value), 1-based, for 4 cells.
write_small_cellranger_directory <- function(directory, entries) {
    dir.create(directory, showWarnings = FALSE)
    writeLines(
        c(
            "gene_1\tGENE1\tGene Expression",
            "grna_1\tgrna_1\tCRISPR Guide Capture",
            "gene_2\tMT-GENE2\tGene Expression",
            "grna_2\tgrna_2\tCRISPR Guide Capture",
            "gene_3\tGENE3\tGene Expression"
        ),
        file.path(directory, "features.tsv")
    )
    writeLines(
        c(
            "%%MatrixMarket matrix coordinate real general",
            paste(5L, 4L, nrow(entries)),
            apply(entries, 1L, paste, collapse = " ")
        ),
        file.path(directory, "matrix.mtx")
    )
    directory
}

small_grna_target_data_frame <- data.frame(
    grna_id = c("grna_1", "grna_2"),
    grna_target = c("target", "non-targeting")
)

import_backed <- function(directories, ...) {
    import_data_from_cellranger(
        directories = directories,
        backed = TRUE,
        directory_to_write = tempfile("backed_import"),
        ...
    ) |>
        suppressMessages()
}

expect_same_import <- function(backed, in_memory) {
    expect_s4_class(get_response_matrix(backed), "DelayedMatrix")
    expect_identical(
        as.matrix(get_response_matrix(backed)),
        as.matrix(get_response_matrix(in_memory))
    )
    expect_identical(
        dimnames(get_response_matrix(backed)),
        dimnames(get_response_matrix(in_memory))
    )
    expect_identical(get_grna_matrix(backed), get_grna_matrix(in_memory))
    expect_identical(
        get_cell_covariates(backed),
        get_cell_covariates(in_memory)
    )
    expect_identical(
        backed@import_grna_assignment_info,
        in_memory@import_grna_assignment_info
    )
    expect_identical(
        backed@grna_target_data_frame,
        in_memory@grna_target_data_frame
    )
}

test_that("a disk-backed Cell Ranger import equals the in-memory import", {
    data(grna_target_data_frame_highmoi)
    directories <- paste0(
        system.file("extdata", package = "sceptre"),
        "/highmoi_example/gem_group_",
        c(1, 2)
    )
    extra_covariates <- data.frame(
        batch = rep(c("x", "y"), length.out = 500L)
    )
    in_memory <- import_data_from_cellranger(
        directories = directories,
        moi = "high",
        grna_target_data_frame = grna_target_data_frame_highmoi,
        extra_covariates = extra_covariates
    ) |>
        suppressMessages()
    backed <- import_backed(
        directories,
        moi = "high",
        grna_target_data_frame = grna_target_data_frame_highmoi,
        extra_covariates = extra_covariates
    )
    expect_same_import(backed, in_memory)
    expect_true(all(
        c("batch", "cellranger_directory") %in%
            colnames(get_cell_covariates(backed))
    ))

    # without barcode files
    copies <- vapply(
        seq_along(directories),
        function(i) {
            copy <- file.path(tempfile("no_barcodes"))
            dir.create(copy)
            file.copy(
                file.path(
                    directories[i],
                    c("features.tsv.gz", "matrix.mtx.gz")
                ),
                copy
            )
            copy
        },
        character(1)
    )
    in_memory <- import_data_from_cellranger(
        directories = copies,
        moi = "high",
        grna_target_data_frame = grna_target_data_frame_highmoi
    ) |>
        suppressMessages()
    backed <- import_backed(
        copies,
        moi = "high",
        grna_target_data_frame = grna_target_data_frame_highmoi
    )
    expect_same_import(backed, in_memory)
})

test_that("a disk-backed import computes p_mito and reads unsorted entries", {
    directory <- write_small_cellranger_directory(
        tempfile("small_cellranger"),
        rbind(
            c(5, 4, 2),
            c(1, 1, 3),
            c(3, 1, 1),
            c(2, 2, 4),
            c(1, 4, 1),
            c(4, 3, 5),
            c(3, 2, 2),
            c(5, 1, 1)
        )
    )
    in_memory <- import_data_from_cellranger(
        directories = directory,
        moi = "low",
        grna_target_data_frame = small_grna_target_data_frame
    ) |>
        suppressMessages()
    backed <- import_backed(
        directory,
        moi = "low",
        grna_target_data_frame = small_grna_target_data_frame
    )
    expect_same_import(backed, in_memory)
    expect_true("response_p_mito" %in% colnames(get_cell_covariates(backed)))
})

test_that("a disk-backed import rejects what the in-memory import rejects", {
    non_integer <- write_small_cellranger_directory(
        tempfile("non_integer"),
        rbind(c(1, 1, 3), c(3, 2, 1.5), c(2, 1, 1), c(4, 2, 2))
    )
    expect_error(
        import_backed(
            non_integer,
            moi = "low",
            grna_target_data_frame = small_grna_target_data_frame
        ),
        "must contain integer counts"
    )
    repeated <- write_small_cellranger_directory(
        tempfile("repeated"),
        rbind(c(1, 1, 3), c(3, 2, 1), c(2, 1, 1), c(4, 2, 2), c(1, 1, 2))
    )
    expect_error(
        import_backed(
            repeated,
            moi = "low",
            grna_target_data_frame = small_grna_target_data_frame
        ),
        "appears more than once"
    )
    expect_error(
        suppressMessages(import_data_from_cellranger(
            directories = repeated,
            moi = "low",
            grna_target_data_frame = small_grna_target_data_frame
        )),
        "appears more than once"
    )
})
