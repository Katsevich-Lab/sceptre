# Shared explicit count fixture; numeric values include stored zeros/fractions.
storage_test_matrix <- function(copies = 1L) {
    mat <- Matrix::sparseMatrix(
        i = rep(c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L), copies),
        j = as.integer(unlist(lapply(seq_len(copies) - 1L, function(k) {
            c(1L, 4L, 8L, 11L, 2L, 5L, 8L, 12L) + 12L * k
        }))),
        x = rep(c(0, 2.5, 3.8, 4.2, 2, 1.2, 0, 8), copies),
        dims = c(2L, 12L * copies), repr = "R"
    )
    rownames(mat) <- c("r1", "r2")
    mat
}
