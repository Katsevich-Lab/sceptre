# run_qc(additional_cells_to_remove =) takes integer indices of the cells (the
# columns of the imported matrices). Anything else is an error rather than a
# silent misreading: a logical vector used to be read as the indices 0 and 1,
# removing only the first cell while the QC metrics reported every TRUE cell as
# removed, and out-of-range indices were ignored.

qc_ready_object <- function() {
    data(lowmoi_example_data, package = "sceptre", envir = environment())
    sceptre_object <- import_data(
        response_matrix = lowmoi_example_data$response_matrix,
        grna_matrix = lowmoi_example_data$grna_matrix,
        grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
        extra_covariates = lowmoi_example_data$extra_covariates,
        moi = "low"
    )
    suppressMessages(
        sceptre_object |>
            set_analysis_parameters(
                discovery_pairs = construct_positive_control_pairs(sceptre_object),
                seed = 4
            ) |>
            assign_grnas(print_progress = FALSE)
    )
}


test_that("integer indices remove exactly those cells, each counted once", {
    sceptre_object <- qc_ready_object()
    covariates <- get_cell_covariates(sceptre_object)
    remove <- which(covariates$batch == covariates$batch[1])
    run <- function(cells) {
        suppressMessages(run_qc(sceptre_object, additional_cells_to_remove = cells))
    }
    by_index <- run(remove)
    expect_false(any(remove %in% by_index@cells_in_use))
    n_removed <- by_index@cell_removal_metrics[["n_cells_rm_user_specified"]]
    expect_identical(n_removed, length(remove))
    # repeated and unordered indices are counted once
    repeated <- run(c(rev(remove), remove[1:3]))
    expect_identical(repeated@cells_in_use, by_index@cells_in_use)
    expect_identical(repeated@cell_removal_metrics, by_index@cell_removal_metrics)
    # a double vector of whole numbers is accepted
    expect_identical(run(as.numeric(remove))@cells_in_use, by_index@cells_in_use)
})


test_that("anything but integer indices of cells is an error", {
    sceptre_object <- qc_ready_object()
    n_cells <- ncol(get_response_matrix(sceptre_object))
    cell_ids <- colnames(get_response_matrix(sceptre_object))
    run <- function(cells) {
        suppressMessages(run_qc(sceptre_object, additional_cells_to_remove = cells))
    }
    expect_error(run(c(1L, n_cells + 1L)), "integer indices")
    expect_error(run(c(0L, 2L)), "integer indices")
    expect_error(run(1.5), "integer indices")
    expect_error(run(c(1L, NA)), "integer indices")
    expect_error(run(rep(TRUE, n_cells)), "logical vectors are not accepted")
    expect_error(run(cell_ids[1:3]), "cell IDs")
    expect_error(run(list(1)), "integer indices")
})
