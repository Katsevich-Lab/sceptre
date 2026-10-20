sce_fixture <- function() {
    skip_if_not_installed("SingleCellExperiment")
    x <- refactor_make_inputs()
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(
            counts = Matrix::Matrix(x$response_matrix, sparse = TRUE)
        ),
        colData = S4Vectors::DataFrame(x$extra_covariates),
        rowData = S4Vectors::DataFrame(symbol = c("MT-A", paste0("gene", 2:8)))
    )
    SingleCellExperiment::altExp(sce, "gRNA") <-
        SummarizedExperiment::SummarizedExperiment(
            assays = list(
                counts = Matrix::Matrix(x$low_grna_matrix, sparse = TRUE)
            )
        )
    object <- import_from_sce(
        sce,
        x$grna_target_data_frame,
        "low",
        covariates = c("batch", "z"),
        response_name_column = "symbol"
    )
    list(sce = sce, object = object, x = x)
}

sce_configure <- function(f) {
    set_analysis_parameters(
        f$object,
        discovery_pairs = f$x$discovery_pairs,
        formula_object = ~ batch + z,
        control_group = "complement",
        resampling_approximation = "no_approximation",
        side = "left"
    )
}

test_that("external assignments work with ordinary matrix import and no SCE", {
    x <- refactor_make_inputs()
    object <- import_data(
        x$response_matrix,
        x$low_grna_matrix,
        x$grna_target_data_frame,
        "low"
    )
    object <- set_analysis_parameters(
        object,
        discovery_pairs = x$discovery_pairs,
        formula_object = ~response_n_umis
    )
    calls <- x$low_grna_matrix >= 40
    assigned <- set_grna_assignments(object, calls, "matrix-only caller v1")
    expect_identical(
        unname(as.matrix(get_grna_assignments(assigned))),
        unname(calls)
    )
    expect_identical(get_cell_covariates(assigned), get_cell_covariates(object))
    expect_identical(
        assigned@spec@grna_assignment_hyperparameters$provenance,
        "matrix-only caller v1"
    )
    expect_error(set_grna_assignments(object, calls[, 160:1], "x"), "cell IDs")
})

test_that("vector imports require canonical vector calls and high MOI retains overlaps", {
    f <- sce_fixture()
    metadata <- f$x$grna_target_data_frame
    metadata$vector_id <- c("va", "va", "vb", "vb", paste0("vn", 1:4))
    o <- import_from_sce(f$sce, metadata, "high")
    o <- set_analysis_parameters(
        o,
        discovery_pairs = f$x$discovery_pairs,
        formula_object = ~response_n_umis
    )
    expect_error(
        set_grna_assignments(o, f$x$low_grna_matrix >= 40, "guide calls"),
        "guide and cell"
    )
    canonical <- as.matrix(get_grna_matrix(o)) >= 40
    canonical[1, 2] <- TRUE
    assigned <- set_grna_assignments(o, canonical, "vector caller v1")
    expect_equal(
        unname(as.matrix(get_grna_assignments(assigned))),
        unname(canonical)
    )
    expect_length(assigned@preparation@assignment_exclusions, 0L)
    expect_equal(
        get_cell_covariates(o)$response_n_umis,
        unname(colSums(f$x$response_matrix))
    )
})

test_that("SCE import preserves raw counts, complete-library covariates and source ownership", {
    f <- sce_fixture()
    before <- serialize(f$sce, NULL)
    expect_identical(
        as.matrix(get_response_matrix(f$object)),
        f$x$response_matrix
    )
    expect_equal(
        get_cell_covariates(f$object)$response_n_umis,
        unname(colSums(f$x$response_matrix))
    )
    expect_equal(
        get_cell_covariates(f$object)$response_p_mito,
        unname(f$x$response_matrix[1, ] / colSums(f$x$response_matrix))
    )
    # A pair subset must not redefine depth, including after QC/serialization.
    o <- sce_configure(f)
    calls <- f$x$low_grna_matrix >= 40
    o <- set_grna_assignments(o, calls, "fixture caller v1")
    o <- run_qc(
        o,
        response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1),
        p_mito_threshold = 1,
        n_nonzero_trt_thresh = 1L,
        n_nonzero_cntrl_thresh = 1L,
        additional_cells_to_remove = c(2L, 17L)
    )
    o <- unserialize(serialize(o, NULL))
    expect_equal(
        get_cell_covariates(o)$response_n_umis,
        unname(colSums(f$x$response_matrix))
    )
    expect_identical(serialize(f$sce, NULL), before)
    expect_identical(calls, f$x$low_grna_matrix >= 40)
})

test_that("SCE boundaries reject misalignment and normalized or missing counts", {
    f <- sce_fixture()
    import <- function(s = f$sce, m = f$x$grna_target_data_frame, ...) {
        import_from_sce(s, m, "low", ...)
    }
    expect_error(
        import(m = f$x$grna_target_data_frame[8:1, ]),
        "match all guide rows"
    )
    expect_error(import(covariates = "missing"), "colData")
    expect_error(import(response_assay = "logcounts"), "raw count assay")
    bad <- f$sce
    child <- SingleCellExperiment::altExp(bad, "gRNA", withDimnames = FALSE)
    child <- child[, rev(seq_len(ncol(child)))]
    SingleCellExperiment::altExp(bad, "gRNA", withDimnames = FALSE) <- child
    # Default SCE retrieval masks the child's stored ordering with parent IDs.
    expect_identical(
        colnames(SingleCellExperiment::altExp(bad, "gRNA")),
        colnames(bad)
    )
    expect_error(import(bad), "cell IDs")
    bad <- f$sce
    physical <- SummarizedExperiment::assay(bad, "counts", withDimnames = FALSE)
    physical <- physical[rev(seq_len(nrow(physical))), , drop = FALSE]
    SummarizedExperiment::assay(bad, "counts", withDimnames = FALSE) <- physical
    expect_error(import(bad), "Stored assay feature IDs")
    bad <- f$sce
    SummarizedExperiment::assay(bad, "counts")[1, 1] <- 0.5
    expect_error(import(bad), "raw counts")
    SummarizedExperiment::assay(bad, "counts")[1, 1] <- Inf
    expect_error(import(bad), "finite")
    bad <- f$sce
    rownames(bad)[2] <- rownames(bad)[1]
    expect_error(import(bad), "unique")
    colnames(bad) <- NULL
    expect_error(import(bad), "IDs must be present")
    if (requireNamespace("DelayedArray", quietly = TRUE)) {
        bad <- f$sce
        SummarizedExperiment::assay(bad, "counts") <-
            DelayedArray::DelayedArray(SummarizedExperiment::assay(
                bad,
                "counts"
            ))
        expect_error(import(bad), "in-memory")
    }
})

test_that("external calls are strict, binary and independent of count thresholding", {
    f <- sce_fixture()
    o <- sce_configure(f)
    calls <- f$x$low_grna_matrix >= 40
    calls[1, 2] <- TRUE # deliberately assign a guide with zero guide counts
    before <- serialize(o, NULL)
    assigned <- set_grna_assignments(o, calls, "other caller v2")
    observed <- as.matrix(get_grna_assignments(assigned))
    expect_equal(unname(observed), unname(calls))
    expect_identical(serialize(o, NULL), before)
    expect_equal(assigned@preparation@assignment_exclusions, 2L)
    expect_identical(assigned@spec@grna_assignment_method, "external")
    expect_error(set_grna_assignments(o, calls[, 160:1], "x"), "cell IDs")
    expect_error(set_grna_assignments(o, calls[8:1, ], "x"), "guide and cell")
    bad <- calls * 1
    bad[1, 1] <- 0.8
    expect_error(set_grna_assignments(o, bad, "x"), "binary")
    bad[1, 1] <- NA
    expect_error(set_grna_assignments(o, bad, "x"), "binary")
    expect_error(set_grna_assignments(o, calls, ""), "provenance")
    expect_error(
        set_grna_assignments(f$object, calls, "x"),
        "set_analysis_parameters"
    )
    empty <- calls
    empty[] <- FALSE
    empty_o <- set_grna_assignments(o, empty, "no calls")
    expect_identical(
        empty_o@preparation@assignment_exclusions,
        seq_len(ncol(empty))
    )
    sparse <- methods::as(Matrix::Matrix(calls, sparse = TRUE), "RsparseMatrix")
    sparse_o <- set_grna_assignments(o, sparse, "sparse caller")
    expect_equal(
        unname(as.matrix(get_grna_assignments(sparse_o))),
        unname(calls)
    )
    pattern <- methods::as(sparse, "nMatrix")
    expect_equal(
        unname(as.matrix(get_grna_assignments(
            set_grna_assignments(o, pattern, "pattern caller")
        ))),
        unname(calls)
    )
    explicit_false <- sparse
    explicit_false[1, 1] <- FALSE
    explicit_false_o <- set_grna_assignments(
        o,
        explicit_false,
        "stored false caller"
    )
    expect_false(as.matrix(get_grna_assignments(explicit_false_o))[1, 1])
})

test_that("external assignment changes invalidate QC/results and stale identities are rejected", {
    f <- sce_fixture()
    calls <- f$x$low_grna_matrix >= 40
    o <- set_grna_assignments(sce_configure(f), calls, "caller v1")
    o <- run_qc(
        o,
        response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1),
        p_mito_threshold = 1,
        n_nonzero_trt_thresh = 1L,
        n_nonzero_cntrl_thresh = 1L
    )
    changed <- set_grna_assignments(o, calls, "caller v2")
    expect_false("cell_qc" %in% names(changed@preparation@stages))
    stale <- o
    stale@preparation@guide_assignments[[1]] <- integer()
    expect_error(get_grna_assignments(stale), "not yet been called")
    expect_true("cell_qc" %in% names(o@preparation@stages))
})

test_that("completed results return as owned pair records without changing the SCE", {
    f <- sce_fixture()
    o <- set_grna_assignments(
        sce_configure(f),
        f$x$low_grna_matrix >= 40,
        "caller v1"
    )
    o <- run_qc(
        o,
        response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1),
        p_mito_threshold = 1,
        n_nonzero_trt_thresh = 1L,
        n_nonzero_cntrl_thresh = 1L
    )
    expect_error(
        export_results_to_sce(o, f$sce),
        "result|run_discovery",
        ignore.case = TRUE
    )
    o <- run_discovery_analysis(o, print_progress = FALSE)
    changed <- set_grna_assignments(
        o,
        f$x$low_grna_matrix >= 40,
        "caller v2 after completed analysis"
    )
    expect_error(
        get_result(changed, "run_discovery_analysis"),
        "result|run_discovery|current|QC"
    )
    before <- serialize(f$sce, NULL)
    result_before <- serialize(get_result(o, "run_discovery_analysis"), NULL)
    exported <- export_results_to_sce(o, f$sce)
    expect_identical(serialize(f$sce, NULL), before)
    expect_identical(
        SummarizedExperiment::assays(exported),
        SummarizedExperiment::assays(f$sce)
    )
    expect_identical(
        SummarizedExperiment::colData(exported),
        SummarizedExperiment::colData(f$sce)
    )
    expect_identical(
        SummarizedExperiment::rowData(exported),
        SummarizedExperiment::rowData(f$sce)
    )
    md <- S4Vectors::metadata(exported)$sceptre
    expect_equal(md$cell_covariates, get_cell_covariates(o))
    expect_identical(md$assignment_method, "external")
    expect_equal(
        md$results$run_discovery_analysis,
        get_result(o, "run_discovery_analysis")
    )
    data.table::set(md$results$run_discovery_analysis, j = "p_value", value = 0)
    expect_identical(
        serialize(get_result(o, "run_discovery_analysis"), NULL),
        result_before
    )
    reopened <- unserialize(serialize(exported, NULL))
    expect_equal(
        S4Vectors::metadata(reopened)$sceptre,
        S4Vectors::metadata(exported)$sceptre
    )
    expect_error(export_results_to_sce(o, exported), "already contains")
    expect_error(export_results_to_sce(o, f$sce[, 160:1]), "exactly match")
    bad <- f$sce
    SummarizedExperiment::assay(bad, "counts")[1, 1] <- 999
    expect_error(export_results_to_sce(o, bad), "exactly match")
})
