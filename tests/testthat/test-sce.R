test_that("as_singlecellexperiment() exports the data in the SCE layout", {
    data(grna_target_data_frame_highmoi)
    data(highmoi_example_data)
    sceptre_object <- import_data(
        response_matrix = highmoi_example_data$response_matrix,
        grna_matrix = highmoi_example_data$grna_matrix,
        grna_target_data_frame = grna_target_data_frame_highmoi,
        moi = "high",
        extra_covariates = highmoi_example_data$extra_covariates
    )
    sce <- as_singlecellexperiment(sceptre_object)
    expect_s4_class(sce, "SingleCellExperiment")
    expect_identical(
        as.matrix(SummarizedExperiment::assay(sce, "counts")),
        as.matrix(get_response_matrix(sceptre_object))
    )
    expect_s4_class(SummarizedExperiment::assay(sce, "counts"), "dgCMatrix")
    expect_identical(
        as.data.frame(SummarizedExperiment::colData(sce), optional = TRUE),
        get_cell_covariates(sceptre_object)
    )
    grna_experiment <- SingleCellExperiment::altExp(sce, "gRNA")
    expect_identical(
        as.matrix(SummarizedExperiment::assay(grna_experiment, "counts")),
        as.matrix(get_grna_matrix(sceptre_object))
    )
    targets <- SummarizedExperiment::rowData(grna_experiment)
    idx <- match(rownames(grna_experiment), grna_target_data_frame_highmoi$grna_id)
    expect_identical(
        targets$grna_target,
        grna_target_data_frame_highmoi$grna_target[idx]
    )
    positions <- SummarizedExperiment::rowRanges(grna_experiment)
    expect_s4_class(positions, "GRangesList")
    positioned <- !is.na(grna_target_data_frame_highmoi$chr[idx])
    expect_identical(unname(lengths(positions)), as.integer(positioned))
    expect_identical(
        GenomicRanges::start(unlist(positions)),
        as.integer(grna_target_data_frame_highmoi$start[idx][positioned])
    )
})

test_that("import_data_from_sce() imports what as_singlecellexperiment() exports", {
    data(grna_target_data_frame_highmoi)
    data(highmoi_example_data)
    sceptre_object <- import_data(
        response_matrix = highmoi_example_data$response_matrix,
        grna_matrix = highmoi_example_data$grna_matrix,
        grna_target_data_frame = grna_target_data_frame_highmoi,
        moi = "high",
        extra_covariates = highmoi_example_data$extra_covariates
    )
    reimported <- import_data_from_sce(
        as_singlecellexperiment(sceptre_object),
        moi = "high",
        extra_covariates = "batch"
    )
    expect_identical(
        get_response_matrix(reimported),
        get_response_matrix(sceptre_object)
    )
    expect_identical(get_grna_matrix(reimported), get_grna_matrix(sceptre_object))
    expect_identical(
        get_cell_covariates(reimported),
        get_cell_covariates(sceptre_object)
    )
    original_targets <- sceptre_object@grna_target_data_frame
    targets <- reimported@grna_target_data_frame
    targets <- targets[match(original_targets$grna_id, targets$grna_id), ]
    rownames(targets) <- NULL
    for (column in c("grna_id", "grna_target", "chr")) {
        expect_identical(targets[[column]], original_targets[[column]])
    }
    for (column in c("start", "end")) {
        expect_equal(targets[[column]], original_targets[[column]])
    }
    expect_identical(
        construct_cis_pairs(reimported),
        construct_cis_pairs(sceptre_object)
    )
})

test_that("import_data_from_sce() reads targets without positions", {
    data(lowmoi_example_data)
    grna_experiment <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = lowmoi_example_data$grna_matrix)
    )
    targets <- lowmoi_example_data$grna_target_data_frame
    SummarizedExperiment::rowData(grna_experiment)$grna_target <-
        targets$grna_target[match(rownames(grna_experiment), targets$grna_id)]
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = lowmoi_example_data$response_matrix),
        colData = lowmoi_example_data$extra_covariates,
        altExps = list(gRNA = grna_experiment)
    )
    from_sce <- import_data_from_sce(sce, moi = "low", extra_covariates = "batch")
    direct <- import_data(
        response_matrix = lowmoi_example_data$response_matrix,
        grna_matrix = lowmoi_example_data$grna_matrix,
        grna_target_data_frame = targets,
        extra_covariates = lowmoi_example_data$extra_covariates,
        moi = "low"
    )
    expect_identical(get_response_matrix(from_sce), get_response_matrix(direct))
    expect_identical(get_cell_covariates(from_sce), get_cell_covariates(direct))
    expected_targets <- as.data.frame(direct@grna_target_data_frame)[
        match(rownames(grna_experiment), targets$grna_id),
        c("grna_id", "grna_target")
    ]
    rownames(expected_targets) <- NULL
    expect_identical(
        from_sce@grna_target_data_frame[c("grna_id", "grna_target")],
        expected_targets
    )
    expect_error(
        import_data_from_sce(sce, moi = "low", grna_altexp = "guides"),
        "`grna_altexp` must name"
    )
    expect_error(
        import_data_from_sce(sce, moi = "low", extra_covariates = "lane"),
        "no column `lane`"
    )
})
