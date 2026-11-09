# `response_names` picks out the mitochondrial genes for `response_p_mito`, so
# it must have one element per response. A longer vector made the C++
# covariate code write past the end of a buffer; a shorter one silently
# computed `response_p_mito` from the wrong genes.

test_that("response_names must have one element per response", {
    data(lowmoi_example_data, package = "sceptre", envir = environment())
    d <- lowmoi_example_data
    n <- nrow(d$response_matrix)
    import <- function(response_names, backed = FALSE) {
        import_data(
            response_matrix = d$response_matrix,
            grna_matrix = d$grna_matrix,
            grna_target_data_frame = d$grna_target_data_frame,
            moi = "low",
            response_names = response_names,
            backed = backed,
            directory_to_write = if (backed) withr::local_tempdir(.local_envir = parent.frame(2))
        )
    }
    too_long <- c(rownames(d$response_matrix), paste0("MT-", seq_len(n)))
    too_short <- c("MT-a", "MT-b", "MT-c")
    for (backed in c(FALSE, TRUE)) {
        expect_error(import(too_long, backed), "one element per response")
        expect_error(import(too_short, backed), "one element per response")
    }

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = d$response_matrix)
    )
    grna <- d$grna_matrix
    rownames(grna) <- rownames(d$grna_matrix)
    SingleCellExperiment::altExp(sce, "grna") <-
        SummarizedExperiment::SummarizedExperiment(assays = list(counts = grna))
    expect_error(
        import_data_from_sce(
            sce,
            moi = "low",
            grna_target_data_frame = d$grna_target_data_frame,
            response_names = too_short
        ),
        "one element per response"
    )
})


test_that("response_p_mito comes from the responses named MT-", {
    data(lowmoi_example_data, package = "sceptre", envir = environment())
    d <- lowmoi_example_data
    response_names <- rownames(d$response_matrix)
    mito <- c(5L, 17L, 60L)
    response_names[mito] <- paste0("MT-", response_names[mito])
    sceptre_object <- import_data(
        response_matrix = d$response_matrix,
        grna_matrix = d$grna_matrix,
        grna_target_data_frame = d$grna_target_data_frame,
        moi = "low",
        response_names = response_names
    )
    counts <- as.matrix(d$response_matrix)
    expected <- colSums(counts[mito, , drop = FALSE]) / colSums(counts)
    expect_equal(
        get_cell_covariates(sceptre_object)$response_p_mito,
        unname(expected)
    )
})


test_that("the C++ covariate code refuses an out-of-range gene index", {
    counts <- Matrix::rsparsematrix(10, 4, 0.5, rand.x = function(n) rpois(n, 3) + 1)
    counts <- methods::as(counts, "CsparseMatrix")
    covariates <- function(mt_gene_idxs) {
        compute_cell_covariates_cpp(
            counts@i, counts@p, counts@x, nrow(counts), ncol(counts),
            mt_gene_idxs, TRUE, FALSE
        )
    }
    expect_no_error(covariates(c(1L, 10L)))
    expect_error(covariates(11L), "out of range")
    expect_error(covariates(0L), "out of range")
})
