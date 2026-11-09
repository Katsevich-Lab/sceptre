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
                barcodes <- read_delimited_file(
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
    expected_barcodes <- read_delimited_file(
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
    barcodes <- read_delimited_file(
        file.path(directory, "barcodes.tsv.gz"),
        header = FALSE,
        colClasses = "character"
    )[[1L]]
    expected_barcodes <- c(
        paste0("b1_", barcodes),
        paste0("b2_", barcodes)
    )

    # the same directory twice
    sceptre_object <- suppressMessages(import_data_from_cellranger(
        directories = c(directory, directory),
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
    barcodes <- read_delimited_file(
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
            read_delimited_file(
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
                barcodes <- read_delimited_file(
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


test_that("Cell Ranger import records the directory of each cell in batch", {
    inputs <- get_cellranger_test_inputs()
    n_cells <- vapply(
        inputs$directories,
        function(directory) {
            nrow(read_delimited_file(
                file.path(directory, "barcodes.tsv.gz"),
                header = FALSE
            ))
        },
        FUN.VALUE = integer(1)
    )

    # two directories: `batch` is b1, b2 in the order of `directories`, and
    # the user's other columns are kept
    user_lane <- rep(c("lane_A", "lane_B", "lane_C"), length.out = sum(n_cells))
    sceptre_object <- suppressMessages(import_data_from_cellranger(
        directories = inputs$directories,
        moi = "high",
        grna_target_data_frame = inputs$grna_target_data_frame,
        extra_covariates = data.frame(lane = user_lane)
    ))
    covariates <- get_cell_covariates(sceptre_object)
    expect_identical(as.character(covariates$lane), user_lane)
    expect_identical(
        covariates$batch,
        factor(rep(c("b1", "b2"), n_cells), levels = c("b1", "b2"))
    )

    # the levels follow the order of `directories`, also beyond nine
    sceptre_object <- suppressMessages(import_data_from_cellranger(
        directories = rep(inputs$directories[1L], 10L),
        moi = "high",
        grna_target_data_frame = inputs$grna_target_data_frame
    ))
    expect_identical(
        levels(get_cell_covariates(sceptre_object)$batch),
        paste0("b", 1:10)
    )

    # with two directories, a user's column named batch is an error, not
    # overwritten
    expect_error(
        import_data_from_cellranger(
            directories = inputs$directories,
            moi = "high",
            grna_target_data_frame = inputs$grna_target_data_frame,
            extra_covariates = data.frame(
                batch = seq_len(sum(n_cells))
            )
        ),
        regexp = "`extra_covariates` contains a column named `batch`",
        fixed = TRUE
    )

    # one directory: the user's batch is kept as supplied
    user_batch <- rep(c("lane_A", "lane_B"), length.out = n_cells[1L])
    sceptre_object <- suppressMessages(import_data_from_cellranger(
        directories = inputs$directories[1L],
        moi = "high",
        grna_target_data_frame = inputs$grna_target_data_frame,
        extra_covariates = data.frame(batch = user_batch)
    ))
    expect_identical(
        as.character(get_cell_covariates(sceptre_object)$batch),
        user_batch
    )
})


test_that("import_data rejects non-integer response counts", {
    data("lowmoi_example_data")
    for (response_matrix in list(
        lowmoi_example_data$response_matrix * 1.5,
        as.matrix(lowmoi_example_data$response_matrix) * 1.5
    )) {
        expect_error(
            import_data(
                response_matrix = response_matrix,
                grna_matrix = lowmoi_example_data$grna_matrix,
                grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
                extra_covariates = lowmoi_example_data$extra_covariates,
                moi = "low"
            ),
            regexp = "The `response_matrix` must contain integer counts",
            fixed = TRUE
        )
    }
})


test_that("import_data rejects missing, negative and infinite counts", {
    data("lowmoi_example_data")
    for (case in list(
        list(value = NA, message = "missing (NA) entries"),
        list(value = -3, message = "negative entries"),
        list(value = Inf, message = "infinite entries")
    )) {
        response_matrix <- lowmoi_example_data$response_matrix
        response_matrix[1L, 1L] <- case$value
        for (input in list(
            response_matrix,
            methods::as(response_matrix, "CsparseMatrix"),
            methods::as(response_matrix, "RsparseMatrix"),
            methods::as(response_matrix, "TsparseMatrix")
        )) {
            expect_error(
                import_data(
                    response_matrix = input,
                    grna_matrix = lowmoi_example_data$grna_matrix,
                    grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
                    moi = "low"
                ),
                regexp = case$message,
                fixed = TRUE
            )
        }
    }
    grna_matrix <- as.matrix(lowmoi_example_data$grna_matrix)
    storage.mode(grna_matrix) <- "integer"
    grna_matrix[1L, 1L] <- -1L
    expect_error(
        import_data(
            response_matrix = lowmoi_example_data$response_matrix,
            grna_matrix = grna_matrix,
            grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
            moi = "low"
        ),
        regexp = "The `grna_matrix` must contain counts",
        fixed = TRUE
    )
    # a logical matrix passes the check unless it has missing entries
    grna_matrix <- as.matrix(lowmoi_example_data$grna_matrix) > 0
    grna_matrix[1L, 1L] <- NA
    expect_error(
        import_data(
            response_matrix = lowmoi_example_data$response_matrix,
            grna_matrix = grna_matrix,
            grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
            moi = "low"
        ),
        regexp = "The `grna_matrix` must contain counts (for example, UMI counts); it has missing (NA) entries",
        fixed = TRUE
    )
})


test_that("import_data rejects non-integer gRNA counts", {
    data("lowmoi_example_data")
    expect_error(
        import_data(
            response_matrix = lowmoi_example_data$response_matrix,
            grna_matrix = lowmoi_example_data$grna_matrix * 1.5,
            grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
            extra_covariates = lowmoi_example_data$extra_covariates,
            moi = "low"
        ),
        regexp = "The `grna_matrix` must contain integer counts",
        fixed = TRUE
    )
    # 0/1 gRNA matrices are accepted
    expect_s4_class(
        suppressMessages(import_data(
            response_matrix = lowmoi_example_data$response_matrix,
            grna_matrix = (lowmoi_example_data$grna_matrix >= 1) * 1,
            grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
            extra_covariates = lowmoi_example_data$extra_covariates,
            moi = "low"
        )),
        "sceptre_object"
    )
})


