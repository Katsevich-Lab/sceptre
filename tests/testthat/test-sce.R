test_that("as_sce() exports the data in the SCE layout", {
    data(grna_target_data_frame_highmoi)
    data(highmoi_example_data)
    sceptre_object <- import_data(
        response_matrix = highmoi_example_data$response_matrix,
        grna_matrix = highmoi_example_data$grna_matrix,
        grna_target_data_frame = grna_target_data_frame_highmoi,
        moi = "high",
        extra_covariates = highmoi_example_data$extra_covariates
    )
    sce <- as_sce(sceptre_object)
    expect_s4_class(sce, "SingleCellExperiment")
    expect_identical(
        as.matrix(SummarizedExperiment::assay(sce, "counts")),
        response_values(sceptre_object)
    )
    expect_s4_class(SummarizedExperiment::assay(sce, "counts"), "dgCMatrix")
    expect_identical(
        as.data.frame(SummarizedExperiment::colData(sce), optional = TRUE),
        get_cell_covariates(sceptre_object)
    )
    grna_experiment <- SingleCellExperiment::altExp(sce, "grna")
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
    # the target regions are rowData columns; rowRanges (the positions of the
    # gRNAs themselves) and the response positions are empty
    for (column in c("chr", "start", "end")) {
        expect_identical(
            targets[[paste0("target_", column)]],
            grna_target_data_frame_highmoi[[column]][idx]
        )
    }
    expect_false(any(c("chr", "start", "end") %in% colnames(targets)))
    expect_true(all(
        lengths(SummarizedExperiment::rowRanges(grna_experiment)) == 0L
    ))
    expect_true(all(lengths(SummarizedExperiment::rowRanges(sce)) == 0L))
    positioned <- !is.na(targets$target_chr)
    target_ranges <- GenomicRanges::makeGRangesFromDataFrame(
        targets[positioned, ],
        keep.extra.columns = TRUE
    )
    expect_identical(
        GenomicRanges::start(target_ranges),
        as.integer(grna_target_data_frame_highmoi$start[idx][positioned])
    )
    expect_identical(
        target_ranges$grna_id,
        rownames(grna_experiment)[positioned]
    )
})

test_that("import_data_from_sce() imports what as_sce() exports", {
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
        as_sce(sceptre_object),
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

test_that("import_data_from_sce() keeps the response positions", {
    data(grna_target_data_frame_highmoi)
    data(highmoi_example_data)
    data(gene_position_data_frame_grch38)
    sce <- as_sce(import_data(
        response_matrix = highmoi_example_data$response_matrix,
        grna_matrix = highmoi_example_data$grna_matrix,
        grna_target_data_frame = grna_target_data_frame_highmoi,
        moi = "high"
    ))
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
    SummarizedExperiment::rowData(sce)$symbol <- highmoi_example_data$gene_names
    sceptre_object <- import_data_from_sce(sce, moi = "high")
    expected_ranges <- gene_ranges
    expect_identical(
        SummarizedExperiment::rowRanges(as_sce(sceptre_object)),
        expected_ranges
    )
    # they survive saving and loading
    directory <- withr::local_tempdir()
    write_sceptre_object(sceptre_object, directory)
    expect_identical(
        SummarizedExperiment::rowRanges(as_sce(read_sceptre_object(directory))),
        expected_ranges
    )
    # they are used in cis pairing when passed, and give the pairs of the
    # packaged positions (in the order of the genes in `sce`)
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
    data(grna_target_data_frame_highmoi)
    data(highmoi_example_data)
    sceptre_object <- import_data(
        response_matrix = highmoi_example_data$response_matrix,
        grna_matrix = highmoi_example_data$grna_matrix,
        grna_target_data_frame = grna_target_data_frame_highmoi,
        moi = "high"
    )
    sce <- as_sce(sceptre_object)
    grna_experiment <- SingleCellExperiment::altExp(sce, "grna")
    targets <- SummarizedExperiment::rowData(grna_experiment)
    colnames(targets) <- sub("^target_", "", colnames(targets))
    SummarizedExperiment::rowData(grna_experiment) <- targets
    SingleCellExperiment::altExp(sce, "grna") <- grna_experiment
    expect_identical(
        construct_cis_pairs(import_data_from_sce(sce, moi = "high")),
        construct_cis_pairs(sceptre_object)
    )
    SummarizedExperiment::rowData(grna_experiment)$target_chr <- targets$chr
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


test_that("as(x, \"SingleCellExperiment\") is as_sce(x)", {
    data(lowmoi_example_data)
    sceptre_object <- suppressMessages(import_data(
        response_matrix = lowmoi_example_data$response_matrix,
        grna_matrix = lowmoi_example_data$grna_matrix,
        grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
        moi = "low"
    ))
    expect_identical(
        methods::as(sceptre_object, "SingleCellExperiment"),
        as_sce(sceptre_object)
    )
})


test_that("as_sce() exports the target regions as given", {
    data(lowmoi_example_data)
    targets <- lowmoi_example_data$grna_target_data_frame
    targets$chr <- ifelse(targets$grna_target == "non-targeting", "", "chr1")
    targets$start <- 100L
    targets$end <- 120L
    sceptre_object <- suppressMessages(import_data(
        response_matrix = lowmoi_example_data$response_matrix,
        grna_matrix = lowmoi_example_data$grna_matrix,
        grna_target_data_frame = targets,
        moi = "low"
    ))
    exported <- SummarizedExperiment::rowData(
        SingleCellExperiment::altExp(as_sce(sceptre_object), "grna")
    )
    idx <- match(exported$grna_id, targets$grna_id)
    expect_identical(exported$target_chr, targets$chr[idx])
    expect_identical(exported$target_start, targets$start[idx])
})
