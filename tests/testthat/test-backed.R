backed_test_inputs <- function() {
    data(lowmoi_example_data)
    list(
        response_matrix = lowmoi_example_data$response_matrix,
        grna_matrix = lowmoi_example_data$grna_matrix,
        grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
        extra_covariates = lowmoi_example_data$extra_covariates,
        moi = "low"
    )
}

run_backed_test_analysis <- function(sceptre_object) {
    positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
    discovery_pairs <- construct_trans_pairs(
        sceptre_object,
        positive_control_pairs,
        "pc_pairs"
    )[seq_len(100L), ]
    sceptre_object |>
        set_analysis_parameters(
            discovery_pairs = discovery_pairs,
            positive_control_pairs = positive_control_pairs,
            resampling_mechanism = "permutations",
            seed = 4L
        ) |>
        assign_grnas() |>
        run_qc() |>
        run_calibration_check(print_progress = FALSE) |>
        run_power_check(print_progress = FALSE) |>
        run_discovery_analysis(print_progress = FALSE) |>
        suppressMessages()
}

expect_same_results <- function(actual, expected) {
    for (analysis in c(
        "run_calibration_check",
        "run_power_check",
        "run_discovery_analysis"
    )) {
        expect_identical(
            as.data.frame(get_result(actual, analysis)),
            as.data.frame(get_result(expected, analysis))
        )
    }
}

test_that("a disk-backed analysis gives the in-memory results", {
    inputs <- backed_test_inputs()
    in_memory <- do.call(import_data, inputs)
    directory <- file.path(tempdir(), "backed analysis")
    backed <- do.call(
        import_data,
        c(inputs, backed = TRUE, directory_to_write = directory)
    )
    expect_s4_class(get_response_matrix(backed), "DelayedMatrix")
    expect_true(file.exists(file.path(directory, "screen.h5mu")))
    expect_identical(
        as.matrix(get_response_matrix(backed)),
        as.matrix(get_response_matrix(in_memory))
    )
    expect_identical(get_grna_matrix(backed), get_grna_matrix(in_memory))
    expect_identical(get_cell_covariates(backed), get_cell_covariates(in_memory))
    expect_false(is.null(h5_counts_location(get_response_matrix(backed))))

    # small blocks, so that the counts are read in many blocks of responses
    old_block_size <- DelayedArray::getAutoBlockSize()
    DelayedArray::setAutoBlockSize(2000) |> suppressMessages()
    on.exit(
        DelayedArray::setAutoBlockSize(old_block_size) |> suppressMessages(),
        add = TRUE
    )
    in_memory <- run_backed_test_analysis(in_memory)
    backed <- run_backed_test_analysis(backed)
    expect_same_results(backed, in_memory)
    expect_identical(backed@M_matrix, in_memory@M_matrix)
    expect_identical(
        backed@n_nonzero_tot_vector,
        in_memory@n_nonzero_tot_vector
    )
    expect_identical(
        backed@discovery_pairs_with_info,
        in_memory@discovery_pairs_with_info
    )
})

test_that("a disk-backed analysis requires permutations", {
    sceptre_object <- do.call(
        import_data,
        c(
            backed_test_inputs(),
            backed = TRUE,
            directory_to_write = file.path(tempdir(), "backed crt")
        )
    )
    expect_error(
        set_analysis_parameters(sceptre_object, resampling_mechanism = "crt"),
        "must be set to 'permutations'"
    )
})

test_that("write_sceptre_object() and read_sceptre_object() round trip", {
    in_memory <- do.call(import_data, backed_test_inputs()) |>
        run_backed_test_analysis()
    directory <- file.path(tempdir(), "written analysis")
    write_sceptre_object(in_memory, directory)
    expect_setequal(
        list.files(directory),
        c("screen.h5mu", "sceptre_object.rds")
    )

    reread <- read_sceptre_object(directory)
    expect_identical(reread@version, in_memory@version)
    expect_identical(
        get_response_matrix(reread),
        get_response_matrix(in_memory)
    )
    expect_identical(get_grna_matrix(reread), get_grna_matrix(in_memory))
    expect_identical(
        get_cell_covariates(reread),
        get_cell_covariates(in_memory)
    )
    expect_same_results(reread, in_memory)

    reread_backed <- read_sceptre_object(directory, backed = TRUE)
    expect_s4_class(get_response_matrix(reread_backed), "DelayedMatrix")
    expect_identical(
        as.matrix(get_response_matrix(reread_backed)),
        as.matrix(get_response_matrix(in_memory))
    )
    expect_same_results(reread_backed, in_memory)
})

test_that("a disk-backed object is written without copying its counts", {
    directory <- file.path(tempdir(), "backed in place")
    backed <- do.call(
        import_data,
        c(backed_test_inputs(), backed = TRUE, directory_to_write = directory)
    )
    screen_file <- file.path(directory, "screen.h5mu")
    modified <- file.mtime(screen_file)
    Sys.sleep(1.1)
    write_sceptre_object(backed, directory)
    expect_identical(file.mtime(screen_file), modified)
    reread <- read_sceptre_object(directory, backed = TRUE)
    expect_identical(
        as.matrix(get_response_matrix(reread)),
        as.matrix(get_response_matrix(backed))
    )

    # written elsewhere, the counts are copied from the backing file
    other_directory <- file.path(tempdir(), "backed copied")
    write_sceptre_object(backed, other_directory)
    in_memory <- do.call(import_data, backed_test_inputs())
    expect_identical(
        get_response_matrix(read_sceptre_object(other_directory)),
        get_response_matrix(in_memory)
    )
})

test_that("read_sceptre_object() checks the cells of screen.h5mu", {
    inputs <- backed_test_inputs()
    original <- suppressMessages(do.call(import_data, inputs))
    order <- rev(seq_len(ncol(inputs$response_matrix)))
    inputs$response_matrix <- inputs$response_matrix[, order]
    inputs$grna_matrix <- inputs$grna_matrix[, order]
    inputs$extra_covariates <- inputs$extra_covariates[order, , drop = FALSE]
    reordered <- suppressMessages(do.call(import_data, inputs))
    directory <- file.path(tempdir(), "original cells")
    other <- file.path(tempdir(), "reordered cells")
    write_sceptre_object(original, directory)
    write_sceptre_object(reordered, other)
    file.copy(
        file.path(other, "screen.h5mu"),
        file.path(directory, "screen.h5mu"),
        overwrite = TRUE
    )
    expect_error(
        read_sceptre_object(directory),
        "does not hold the data"
    )
})

test_that("read_sceptre_object() checks the directory", {
    directory <- file.path(tempdir(), "incomplete analysis")
    dir.create(directory, showWarnings = FALSE)
    expect_error(read_sceptre_object(directory), "must contain the files")
    expect_error(read_sceptre_object(directory, backed = NA), "`backed`")
})

test_that("use_ondisc and the ondisc functions are deprecated", {
    directory <- file.path(tempdir(), "deprecated ondisc")
    expect_warning(
        sceptre_object <- do.call(
            import_data,
            c(
                backed_test_inputs(),
                use_ondisc = TRUE,
                directory_to_write = directory
            )
        ),
        "`use_ondisc` is deprecated"
    )
    expect_s4_class(get_response_matrix(sceptre_object), "DelayedMatrix")

    expect_warning(
        write_ondisc_backed_sceptre_object(sceptre_object, directory),
        "deprecated"
    )
    expect_warning(
        reread <- read_ondisc_backed_sceptre_object(
            file.path(directory, "sceptre_object.rds")
        ),
        "deprecated"
    )
    expect_s4_class(get_response_matrix(reread), "DelayedMatrix")
    warnings <- testthat::capture_warnings(
        read_ondisc_backed_sceptre_object(
            file.path(directory, "sceptre_object.rds"),
            response_odm_file_fp = file.path(directory, "response.odm"),
            grna_odm_file_fp = file.path(directory, "grna.odm")
        )
    )
    expect_true(any(grepl("are ignored", warnings)))

    # the files of an ondisc-backed object of an earlier version
    old_directory <- file.path(tempdir(), "old ondisc")
    dir.create(old_directory, showWarnings = FALSE)
    saveRDS(list(), file.path(old_directory, "sceptre_object.rds"))
    expect_error(
        suppressWarnings(read_ondisc_backed_sceptre_object(
            file.path(old_directory, "sceptre_object.rds")
        )),
        "re-import the data"
    )
})

test_that("screen.h5mu requires unique cell IDs", {
    screen <- new_screen(
        response_matrix = matrix(
            1,
            2L,
            2L,
            dimnames = list(c("r1", "r2"), c("c", "c"))
        ),
        grna_matrix = matrix(1, 1L, 2L, dimnames = list("g", c("c", "c")))
    )
    expect_error(
        write_screen(screen, file.path(tempdir(), "duplicated.h5mu")),
        "cell IDs must be unique"
    )
})

test_that("pairwise QC handles an analysis without pairs", {
    inputs <- backed_test_inputs()
    for (backed in c(FALSE, TRUE)) {
        sceptre_object <- do.call(
            import_data,
            c(
                inputs,
                backed = backed,
                directory_to_write = file.path(tempdir(), "backed no pairs")
            )
        ) |>
            set_analysis_parameters(resampling_mechanism = "permutations") |>
            assign_grnas() |>
            run_qc() |>
            suppressMessages()
        expect_identical(nrow(sceptre_object@discovery_pairs_with_info), 0L)
        expect_length(
            sceptre_object@n_nonzero_tot_vector,
            nrow(inputs$response_matrix)
        )
    }
})
