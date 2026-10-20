seurat_fixture <- function(legacy = FALSE) {
    skip_if_not_installed("SeuratObject", "5.0.0")
    x <- refactor_make_inputs()
    response <- methods::as(
        Matrix::Matrix(x$response_matrix, sparse = TRUE),
        "generalMatrix"
    )
    grna <- methods::as(
        Matrix::Matrix(x$low_grna_matrix, sparse = TRUE),
        "generalMatrix"
    )
    make_assay <- if (legacy) {
        SeuratObject::CreateAssayObject
    } else {
        SeuratObject::CreateAssay5Object
    }
    object <- SeuratObject::CreateSeuratObject(
        counts = make_assay(counts = response, min.cells = 0, min.features = 0),
        meta.data = x$extra_covariates
    )
    object[["gRNA"]] <- make_assay(
        counts = grna,
        min.cells = 0,
        min.features = 0
    )
    list(seurat = object, x = x)
}

test_that("Seurat imports complete raw v5 and legacy assays without source mutation", {
    for (legacy in c(FALSE, TRUE)) {
        f <- seurat_fixture(legacy)
        before <- serialize(f$seurat, NULL)
        o <- import_from_seurat(
            f$seurat,
            f$x$grna_target_data_frame,
            "low",
            covariates = c("batch", "z")
        )
        expect_identical(as.matrix(get_response_matrix(o)), f$x$response_matrix)
        expect_identical(as.matrix(get_grna_matrix(o)), f$x$low_grna_matrix)
        expect_equal(
            get_cell_covariates(o)$response_n_umis,
            unname(colSums(f$x$response_matrix))
        )
        expect_identical(serialize(f$seurat, NULL), before)
        s <- f$seurat
        s[["RNA"]][["symbol"]] <- stats::setNames(
            c("MT-A", paste0("gene", 2:8)),
            rownames(s[["RNA"]])
        )
        named <- import_from_seurat(
            s,
            f$x$grna_target_data_frame,
            "low",
            response_name_column = "symbol"
        )
        expect_equal(
            get_cell_covariates(named)$response_p_mito,
            unname(f$x$response_matrix[1, ] / colSums(f$x$response_matrix))
        )
        expect_identical(
            get_cell_covariates(o)$batch,
            f$x$extra_covariates$batch
        )
    }
})

test_that("Seurat rejects split counts and cell-specific or feature-specific subsets", {
    f <- seurat_fixture()
    import <- function(s) {
        import_from_seurat(s, f$x$grna_target_data_frame, "low")
    }
    bad <- f$seurat
    count <- SeuratObject::LayerData(bad[["RNA"]], layer = "counts")
    bad[["RNA"]] <- SeuratObject::CreateAssay5Object(
        counts = list(a = count[, seq_len(80)], b = count[, 81:160]),
        min.cells = 0,
        min.features = 0
    )
    expect_error(import(bad), "unsplit counts")
    bad <- f$seurat
    bad[["gRNA"]] <- subset(bad[["gRNA"]], cells = colnames(bad)[-1L])
    expect_error(import(bad), "all assay features and object cells")
    bad <- f$seurat
    SeuratObject::LayerData(bad[["RNA"]], layer = "data") <- count
    SeuratObject::LayerData(bad[["RNA"]], layer = "counts") <- count[-1L, ]
    expect_error(import(bad), "all assay features and object cells")
    bad <- f$seurat
    bad[["RNA"]] <- SeuratObject::CreateAssay5Object(
        data = count,
        min.cells = 0,
        min.features = 0
    )
    expect_error(import(bad), "unsplit counts")
})

test_that("Seurat does not infer guide metadata or mask named physical permutations", {
    f <- seurat_fixture()
    expect_error(
        import_from_seurat(f$seurat, f$x$grna_target_data_frame[8:1, ], "low"),
        "match all guide rows"
    )
    bad <- f$seurat
    # Deliberate corruption tests a valid named layer with stale feature maps.
    assay <- bad[["RNA"]]
    raw <- SeuratObject::LayerData(assay, layer = "counts", fast = TRUE)
    dimnames(raw) <- list(rev(rownames(assay)), colnames(assay))
    assay@layers$counts <- raw
    bad[["RNA"]] <- assay
    expect_error(
        import_from_seurat(bad, f$x$grna_target_data_frame, "low"),
        "Stored count layer IDs"
    )
    bad <- f$seurat
    raw <- SeuratObject::LayerData(bad[["RNA"]], layer = "counts")
    raw[1, 1] <- 0.25
    SeuratObject::LayerData(bad[["RNA"]], layer = "counts") <- raw
    expect_error(
        import_from_seurat(bad, f$x$grna_target_data_frame, "low"),
        "raw counts"
    )
})

test_that("Seurat returns owned pair results and preserves original metadata and assays", {
    f <- seurat_fixture()
    embedding <- cbind(
        PC_1 = sin(seq_len(ncol(f$seurat))),
        PC_2 = cos(seq_len(ncol(f$seurat)))
    )
    rownames(embedding) <- colnames(f$seurat)
    f$seurat[["pca"]] <- SeuratObject::CreateDimReducObject(
        embeddings = embedding,
        key = "PC_",
        assay = "RNA"
    )
    SeuratObject::Misc(f$seurat, slot = "user_note") <- list(text = "preserve")
    o <- import_from_seurat(f$seurat, f$x$grna_target_data_frame, "low")
    o <- set_analysis_parameters(
        o,
        discovery_pairs = f$x$discovery_pairs,
        formula_object = ~response_n_umis,
        control_group = "complement",
        side = "left",
        resampling_approximation = "no_approximation"
    )
    o <- set_grna_assignments(
        o,
        f$x$low_grna_matrix >= 40,
        "Seurat example caller"
    )
    o <- run_qc(
        o,
        response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1),
        p_mito_threshold = 1,
        n_nonzero_trt_thresh = 1L,
        n_nonzero_cntrl_thresh = 1L
    )
    expect_error(export_results_to_seurat(o, f$seurat), "result")
    o <- run_discovery_analysis(o, print_progress = FALSE)
    before <- serialize(f$seurat, NULL)
    returned <- export_results_to_seurat(o, f$seurat)
    expect_identical(serialize(f$seurat, NULL), before)
    expect_identical(returned[[]], f$seurat[[]])
    expect_identical(returned[["RNA"]], f$seurat[["RNA"]])
    expect_identical(returned[["gRNA"]], f$seurat[["gRNA"]])
    expect_identical(returned[["pca"]], f$seurat[["pca"]])
    expect_identical(
        SeuratObject::Misc(returned, slot = "user_note"),
        SeuratObject::Misc(f$seurat, slot = "user_note")
    )
    md <- SeuratObject::Misc(returned, slot = "sceptre")
    expect_equal(
        md$results$run_discovery_analysis,
        get_result(o, "run_discovery_analysis")
    )
    data.table::set(md$results$run_discovery_analysis, j = "p_value", value = 0)
    expect_false(all(
        get_result(o, "run_discovery_analysis")$p_value == 0,
        na.rm = TRUE
    ))
    expect_error(export_results_to_seurat(o, returned), "already contains")
    bad <- f$seurat
    raw <- SeuratObject::LayerData(bad[["RNA"]], layer = "counts")
    raw[1, 1] <- 999
    SeuratObject::LayerData(bad[["RNA"]], layer = "counts") <- raw
    expect_error(export_results_to_seurat(o, bad), "match the imported library")
    expect_equal(
        SeuratObject::Misc(unserialize(serialize(returned, NULL))),
        SeuratObject::Misc(returned)
    )
})
