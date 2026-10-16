# Test-only original implementation from 21f9ea098b69c444a884a49bc254de711968e3b5
# R/import_functs.R; function bytes below are unchanged. Never used in production.
collapse_grna_matrix <- function(grna_matrix, grna_target_data_frame) {
    # temporary function; we may rewrite this function to improve its
    # speed and memory efficiency if multiguide vector data become more
    # prevalent.
    vector_ids <- grna_target_data_frame$vector_id |> unique()
    grna_matrix <- as.matrix(grna_matrix)
    n_cells <- ncol(grna_matrix)
    grna_matrix_collapsed <- vapply(
        vector_ids,
        function(curr_vector_id) {
            curr_grna_ids <- grna_target_data_frame |>
                dplyr::filter(vector_id == curr_vector_id) |>
                dplyr::pull(grna_id)
            Matrix::colSums(grna_matrix[curr_grna_ids, , drop = FALSE])
        },
        FUN.VALUE = numeric(n_cells)
    ) |>
        t()
    return(grna_matrix_collapsed)
}
