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
        Map(
            function(directory, batch_id) {
                barcodes <- data.table::fread(
                    file.path(directory, "barcodes.tsv.gz"),
                    header = FALSE,
                    colClasses = "character"
                )[[1L]]
                paste0("b", batch_id, "_", barcodes)
            },
            inputs$directories,
            seq_along(inputs$directories)
        ),
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


test_that("single-directory Cell Ranger imports retain raw cell barcodes", {
    inputs <- get_cellranger_test_inputs()
    directory <- inputs$directories[1L]
    expected_barcodes <- data.table::fread(
        file.path(directory, "barcodes.tsv.gz"),
        header = FALSE,
        colClasses = "character"
    )[[1L]]

    sceptre_object <- suppressMessages(import_data_from_cellranger(
        directories = directory,
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


test_that("batch prefixes disambiguate duplicate Cell Ranger barcodes", {
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


test_that("multi-directory barcodes agree with prefixed extra covariates", {
    inputs <- get_cellranger_test_inputs()
    raw_barcodes <- unlist(
        lapply(inputs$directories, function(directory) {
            data.table::fread(
                file.path(directory, "barcodes.tsv.gz"),
                header = FALSE,
                colClasses = "character"
            )[[1L]]
        }),
        use.names = FALSE
    )
    prefixed_barcodes <- unlist(
        Map(
            function(directory, batch_id) {
                barcodes <- data.table::fread(
                    file.path(directory, "barcodes.tsv.gz"),
                    header = FALSE,
                    colClasses = "character"
                )[[1L]]
                paste0("b", batch_id, "_", barcodes)
            },
            inputs$directories,
            seq_along(inputs$directories)
        ),
        use.names = FALSE
    )
    extra_covariates <- data.frame(
        covariate = seq_along(prefixed_barcodes),
        row.names = prefixed_barcodes
    )

    sceptre_object <- suppressMessages(import_data_from_cellranger(
        directories = inputs$directories,
        moi = "high",
        grna_target_data_frame = inputs$grna_target_data_frame,
        extra_covariates = extra_covariates
    ))
    expect_identical(
        rownames(get_cell_covariates(sceptre_object)),
        prefixed_barcodes
    )

    rownames(extra_covariates) <- raw_barcodes
    expect_error(
        suppressMessages(import_data_from_cellranger(
            directories = inputs$directories,
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

test_that("Cell Ranger import gives each vector one row spanning its gRNAs", {
    inputs <- get_cellranger_test_inputs()
    targets <- inputs$grna_target_data_frame
    # one vector per target, holding its gRNAs, whose positions differ
    targets$vector_id <- ifelse(
        targets$grna_target == "non-targeting",
        targets$grna_id,
        paste0(targets$grna_target, "_v1")
    )
    repeated <- rbind(targets, targets[1L, ])
    modes <- if (requireNamespace("ondisc", quietly = TRUE)) {
        c(FALSE, TRUE)
    } else {
        FALSE
    }
    for (use_ondisc in modes) {
        import <- function(grna_target_data_frame) {
            suppressMessages(import_data_from_cellranger(
                directories = inputs$directories,
                moi = "high",
                grna_target_data_frame = grna_target_data_frame,
                use_ondisc = use_ondisc,
                directory_to_write = if (use_ondisc) {
                    withr::local_tempdir(.local_envir = parent.frame(2L))
                }
            ))
        }
        object <- import(targets)
        vectors <- object@grna_target_data_frame
        expect_identical(vectors$grna_id, unique(targets$vector_id))
        expect_identical(
            vectors$start,
            as.vector(tapply(targets$start, targets$vector_id, min)[
                vectors$grna_id
            ])
        )
        expect_identical(
            vectors$end,
            as.vector(tapply(targets$end, targets$vector_id, max)[
                vectors$grna_id
            ])
        )
        # with one row per vector, singleton integration tests each pair once
        object <- set_analysis_parameters(
            object,
            discovery_pairs = data.frame(
                grna_target = vectors$grna_target[1:2],
                response_id = rownames(get_response_matrix(object))[1:2]
            ),
            grna_integration_strategy = "singleton",
            resampling_mechanism = "permutations"
        )
        expect_identical(anyDuplicated(object@discovery_pairs_with_info), 0L)
        expect_error(
            import(repeated),
            paste("These appear more than once:", targets$grna_id[1L]),
            fixed = TRUE
        )
    }
})

test_that("a vector's gRNAs must agree in every column but start and end", {
    targets <- data.frame(
        grna_id = c("g1", "g2", "g3"),
        grna_target = c("A", "A", "B"),
        chr = c("chr1", "chr2", "chr1"),
        start = c(10, 50, 90),
        end = c(30, 70, 110),
        vector_id = c("v1", "v1", "v2")
    )
    expect_error(
        collapse_grna_target_data_frame(targets),
        "The gRNAs of vector `v1` must agree in `chr`",
        fixed = TRUE
    )
    targets$chr <- "chr1"
    expect_identical(
        collapse_grna_target_data_frame(targets),
        data.frame(
            grna_target = c("A", "B"),
            chr = "chr1",
            start = c(10, 90),
            end = c(70, 110),
            grna_id = c("v1", "v2")
        )
    )
})
