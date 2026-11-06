# The high-MOI example data as a SingleCellExperiment in the layout that
# import_data_from_sce() reads, the target regions in the rowData columns
# `target_chr`, `target_start` and `target_end`, or, with `prefix = ""`, in
# `chr`, `start` and `end`
highmoi_sce <- function(prefix = "target_") {
    data(grna_target_data_frame_highmoi, envir = environment())
    data(highmoi_example_data, envir = environment())
    grna_experiment <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = highmoi_example_data$grna_matrix)
    )
    targets <- grna_target_data_frame_highmoi[match(
        rownames(grna_experiment),
        grna_target_data_frame_highmoi$grna_id
    ), ]
    row_data <- SummarizedExperiment::rowData(grna_experiment)
    row_data$grna_target <- targets$grna_target
    for (column in c("chr", "start", "end")) {
        row_data[[paste0(prefix, column)]] <- targets[[column]]
    }
    SummarizedExperiment::rowData(grna_experiment) <- row_data
    SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = highmoi_example_data$response_matrix),
        colData = highmoi_example_data$extra_covariates,
        altExps = list(grna = grna_experiment)
    )
}

test_that("import_data_from_sce() imports the data of the SCE layout", {
    data(grna_target_data_frame_highmoi)
    data(highmoi_example_data)
    direct <- import_data(
        response_matrix = highmoi_example_data$response_matrix,
        grna_matrix = highmoi_example_data$grna_matrix,
        grna_target_data_frame = grna_target_data_frame_highmoi,
        moi = "high",
        extra_covariates = highmoi_example_data$extra_covariates
    )
    from_sce <- import_data_from_sce(
        highmoi_sce(),
        moi = "high",
        extra_covariates = "batch"
    )
    expect_identical(
        get_response_matrix(from_sce),
        get_response_matrix(direct)
    )
    expect_identical(get_grna_matrix(from_sce), get_grna_matrix(direct))
    expect_identical(
        get_cell_covariates(from_sce),
        get_cell_covariates(direct)
    )
    # the target region columns become the columns chr, start and end
    original_targets <- direct@grna_target_data_frame
    targets <- from_sce@grna_target_data_frame
    targets <- targets[match(original_targets$grna_id, targets$grna_id), ]
    rownames(targets) <- NULL
    for (column in c("grna_id", "grna_target", "chr")) {
        expect_identical(targets[[column]], original_targets[[column]])
    }
    for (column in c("start", "end")) {
        expect_equal(targets[[column]], original_targets[[column]])
    }
    expect_identical(construct_cis_pairs(from_sce), construct_cis_pairs(direct))
})

test_that("import_data_from_sce() does not use the response positions", {
    data(gene_position_data_frame_grch38)
    sce <- highmoi_sce()
    without_positions <- import_data_from_sce(sce, moi = "high")
    # the genes at their packaged positions, with a rowData column
    positions <- gene_position_data_frame_grch38[
        match(rownames(sce), gene_position_data_frame_grch38$response_id),
    ]
    gene_ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
        chr = positions$chr,
        start = positions$position,
        end = positions$position
    ))
    names(gene_ranges) <- rownames(sce)
    SummarizedExperiment::rowRanges(sce) <- gene_ranges
    SummarizedExperiment::rowData(sce)$symbol <- rownames(sce)
    sceptre_object <- import_data_from_sce(sce, moi = "high")
    expect_identical(
        get_response_matrix(sceptre_object),
        get_response_matrix(without_positions)
    )
    expect_identical(
        construct_cis_pairs(sceptre_object, distance_threshold = 5e6),
        construct_cis_pairs(without_positions, distance_threshold = 5e6)
    )
    # passed to construct_cis_pairs(), they give the pairs of the packaged
    # positions (in the order of the genes in `sce`)
    sorted_pairs <- function(pairs) {
        pairs <- pairs[order(pairs$grna_target, pairs$response_id), ]
        rownames(pairs) <- NULL
        pairs
    }
    expected_pairs <- construct_cis_pairs(
        sceptre_object,
        distance_threshold = 5e6
    )
    expect_gt(nrow(expected_pairs), 0L)
    expect_identical(
        sorted_pairs(construct_cis_pairs(
            sceptre_object,
            distance_threshold = 5e6,
            response_position_data_frame = SummarizedExperiment::rowRanges(sce)
        )),
        sorted_pairs(expected_pairs)
    )
})

test_that("import_data_from_sce() reads target regions named chr, start, end", {
    sce <- highmoi_sce(prefix = "")
    expect_identical(
        construct_cis_pairs(import_data_from_sce(sce, moi = "high")),
        construct_cis_pairs(import_data_from_sce(highmoi_sce(), moi = "high"))
    )
    grna_experiment <- SingleCellExperiment::altExp(sce, "grna")
    SummarizedExperiment::rowData(grna_experiment)$target_chr <-
        SummarizedExperiment::rowData(grna_experiment)$chr
    SingleCellExperiment::altExp(sce, "grna") <- grna_experiment
    expect_error(
        import_data_from_sce(sce, moi = "high"),
        "not both",
        fixed = TRUE
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
        altExps = list(grna = grna_experiment)
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
    backed <- import_data_from_sce(
        sce,
        moi = "low",
        extra_covariates = "batch",
        backed = TRUE,
        directory_to_write = file.path(tempdir(), "sce backed")
    )
    expect_s4_class(get_response_matrix(backed), "DelayedMatrix")
    expect_identical(
        response_values(backed),
        response_values(direct)
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

