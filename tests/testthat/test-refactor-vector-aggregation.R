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
        # Valid triplets permit duplicate coordinates; their Matrix value
        # aggregation differs from repeated metadata selections. Stored FALSE
        # is zero here, regardless of storage's separate position-count contract.
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

test_that("vector sums retain current large-count and fractional rounding behavior", {
    counts <- matrix(c(2^53, 1, 1, 0, .1, .2, .3, 0), 4L,
                     dimnames = list(paste0("g", 1:4), c("large", "fractional")))
    metadata <- data.frame(grna_id = paste0("g", 1:3), vector_id = "v", grna_target = "t")
    # Verified original Matrix::colSums double accumulation: unit increments at
    # 2^53 are lost, and .1+.2+.3 rounds above .6. Do not silently compensate it.
    expected <- matrix(c(2^53, .1 + .2 + .3), 1L, dimnames = list("v", colnames(counts)))
    for (input in list(counts, methods::as(Matrix::Matrix(counts, sparse = TRUE), "RsparseMatrix"),
                      methods::as(Matrix::Matrix(counts, sparse = TRUE), "TsparseMatrix"))) {
        expect_identical(refactor_private("collapse_grna_matrix")(input, metadata), expected)
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
            expect_identical(object@grna_target_data_frame,
                data.frame(grna_target = c("targetZ", "targetA"), grna_id = c("vZ", "vA")))
            expect_identical(serialize(list(counts, response, fixture$metadata), NULL), before)
        }
    }
})
