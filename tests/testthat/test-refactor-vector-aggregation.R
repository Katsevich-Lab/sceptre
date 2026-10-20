aggregation_example <- function(logical = FALSE) {
    i <- c(0L, 0L, 1L, 2L, 3L, 1L, 3L)
    j <- c(0L, 0L, 1L, 2L, 0L, 0L, 2L)
    dims <- c(4L, 3L)
    labels <- list(paste0("g", 1:4), c("c-z", "c-a", "c-m"))
    if (logical) {
        counts <- methods::new("lgTMatrix", i = i, j = j,
            x = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
            Dim = dims, Dimnames = labels)
        expected <- rbind(vZ = c(0, 0, 1), vA = c(3, 0, 0))
    } else {
        counts <- methods::new("dgTMatrix", i = i, j = j,
            x = c(.25, .75, 0, 1.5, 2, 3, .5), Dim = dims, Dimnames = labels)
        expected <- rbind(vZ = c(2, 0, 2), vA = c(5, 0, 0))
    }
    colnames(expected) <- labels[[2L]]
    metadata <- data.frame(grna_id = c("g3", "g1", "g1", "g2", "g4"),
        grna_target = c("targetZ", "targetA", "targetA", "targetA", "targetZ"),
        vector_id = c("vZ", "vA", "vA", "vA", "vZ"))
    list(counts = counts, metadata = metadata, expected = expected)
}

test_that("vector aggregation preserves dense output and both duplicate meanings", {
    collapse <- refactor_private("collapse_grna_matrix")
    for (logical in c(FALSE, TRUE)) {
        fixture <- aggregation_example(logical)
        # Matrix can combine duplicate sparse triplet entries. That differs
        # from selecting a guide twice through repeated metadata rows. Here
        # stored FALSE contributes zero; QC's stored-entry counts use a different rule.
        representations <- list(dense = as.matrix(fixture$counts),
            C = methods::as(fixture$counts, "CsparseMatrix"),
            R = methods::as(fixture$counts, "RsparseMatrix"), T = fixture$counts)
        for (counts in representations) {
            if (methods::is(counts, "sparseMatrix")) expect_true(methods::validObject(counts))
            before <- serialize(list(counts, fixture$metadata), NULL)
            actual <- collapse(counts, fixture$metadata)
            expect_identical(actual, fixture$expected)
            expect_identical(class(actual), c("matrix", "array"))
            expect_type(actual, "double")
            expect_identical(serialize(list(counts, fixture$metadata), NULL), before)
        }
    }
})

# Calculate expected sums directly from dense matrices, keeping repeated
# metadata rows. Use base colSums() to match this platform's numeric rounding.
aggregation_dense_spec <- function(counts, metadata) {
    dense <- as.matrix(counts)
    groups <- unique(metadata$vector_id)
    labels <- if (is.character(groups)) groups else NULL
    out <- matrix(0, length(groups), ncol(dense),
                  dimnames = list(labels, if (length(groups)) colnames(dense) else NULL))
    for (i in seq_along(groups)) {
        # which() drops missing group comparisons. Factor row selectors retain
        # R matrix indexing semantics; character selectors use row labels.
        selected <- which(metadata$vector_id == groups[[i]])
        rows <- metadata$grna_id[selected]
        out[i, ] <- base::colSums(dense[rows, , drop = FALSE])
    }
    # The public implementation retains its legacy one-cell simplification:
    # a 1-by-group matrix, with group names on columns instead of rows.
    if (ncol(dense) == 1L)
        return(matrix(as.numeric(out), 1L, length(groups), dimnames = list(NULL, labels)))
    out
}

test_that("vector sums preserve the platform's dense column-sum contract", {
    counts <- matrix(c(2^53, 1, 1, 0, .1, .2, .3, 0), 4L,
                     dimnames = list(paste0("g", 1:4), c("large", "fractional")))
    metadata <- data.frame(grna_id = paste0("g", 1:3), vector_id = "v", grna_target = "t")
    for (input in list(counts, methods::as(Matrix::Matrix(counts, sparse = TRUE), "RsparseMatrix"),
                      methods::as(Matrix::Matrix(counts, sparse = TRUE), "TsparseMatrix"))) {
        expect_identical(refactor_private("collapse_grna_matrix")(input, metadata),
                         aggregation_dense_spec(input, metadata))
    }
})

test_that("column blocks keep every selected row in the original sum", {
    # 1025 x1025 crosses the 2^20-element target and leaves a final short block.
    n <- 1025L
    counts <- Matrix::sparseMatrix(i = rep(1:3, n), j = rep(seq_len(n), each = 3L),
        x = rep(c(2^53, 1, 1), n), dims = c(n, n),
        dimnames = list(paste0("g", seq_len(n)), paste0("c", seq_len(n))))
    metadata <- data.frame(grna_id = rownames(counts), vector_id = "v", grna_target = "t")
    before <- serialize(counts, NULL)
    expect_identical(refactor_private("collapse_grna_matrix")(counts, metadata),
                     aggregation_dense_spec(counts, metadata))
    expect_identical(serialize(counts, NULL), before)
    for (columns in list(integer(), 1L)) {
        input <- counts[, columns, drop = FALSE]
        expect_identical(refactor_private("collapse_grna_matrix")(input, metadata),
                         aggregation_dense_spec(input, metadata))
    }
})

test_that("vector aggregation agrees with dense algebra on fixed and seeded duplicates", {
    path <- testthat::test_path("fixtures", "aggregation-inputs-6a43fa5.rds")
    stopifnot(identical(digest::digest(file = path, algo = "sha256"),
        "f53c93eabbb30a53d5fed01788ee8ef0ef4447b94ed308540410e21baa624e2d"))
    fixtures <- readRDS(path)
    expect_length(fixtures, 28L)
    withr::with_seed(94201L, {
        for (k in seq_len(200L)) {
            counts <- methods::new("dgTMatrix", i = sample.int(4L, 30L, TRUE) - 1L,
                j = sample.int(3L, 30L, TRUE) - 1L,
                x = sample(c(0, .1, .2, .3, 1, 2, 2^53, 1e-20), 30L, TRUE),
                Dim = c(4L, 3L), Dimnames = list(paste0("g", 1:4), c("z", "a", "m")))
            metadata <- data.frame(grna_id = paste0("g", sample.int(4L, 8L, TRUE)),
                grna_target = "t", vector_id = sample(c("vZ", "vA"), 8L, TRUE))
            fixtures[[paste0("rounding-", k)]] <- list(m = counts, d = metadata)
        }
    })
    for (fixture in fixtures) {
        before <- serialize(fixture, NULL)
        expect_identical(refactor_private("collapse_grna_matrix")(fixture$m, fixture$d),
                         aggregation_dense_spec(fixture$m, fixture$d))
        expect_identical(serialize(fixture, NULL), before)
    }
})

test_that("public in-memory vector imports retain cell order and metadata multiplicity", {
    for (logical in c(FALSE, TRUE)) {
        fixture <- aggregation_example(logical)
        for (representation in c("CsparseMatrix", "RsparseMatrix", "TsparseMatrix")) {
            counts <- methods::as(fixture$counts, representation)
            response <- matrix(c(1, 2, 3, 4, 5, 6), 2L,
                dimnames = list(c("r-z", "r-a"), colnames(counts)))
            before <- serialize(list(counts, response, fixture$metadata), NULL)
            object <- sceptre::import_data(response, counts, fixture$metadata, moi = "high")
            expect_identical(as.matrix(sceptre::get_grna_matrix(object)), fixture$expected)
            expect_s4_class(sceptre::get_grna_matrix(object), "dgRMatrix")
            expect_identical(rownames(sceptre::get_cell_covariates(object)), colnames(counts))
            expect_identical(object@data@grna_target_data_frame,
                data.frame(grna_target = c("targetZ", "targetA"), grna_id = c("vZ", "vA")))
            expect_identical(serialize(list(counts, response, fixture$metadata), NULL), before)
        }
    }
})
