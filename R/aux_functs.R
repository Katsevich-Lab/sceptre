#' Combine perturbations
#'
#'
#'
#' @param perturbation_matrix a perturbation matrix (i.e., a binary matrix of gRNA-to-cell assignments) stored as a sparse matrix (as implemented by the Matrix package) or a dense matrix (as implemented by base R); the row names should be the gRNA IDs
#' @param gRNA_groups_table a data frame with columns `gRNA_id` and `gRNA_group`, mapping each gRNA to its gRNA group
#'
#' @return a "combined" perturbation matrix, where gRNAs within the same gRNA group have been collapsed into a single row
#' @export
#'
#' @examples
#' library(Matrix)
#' data("gRNA_matrix")
#' data("gRNA_groups_table")
#' perturbation_matrix <- threshold_gRNA_matrix(gRNA_matrix)
#' combined_perturbation_matrix <- combine_perturbations(perturbation_matrix, gRNA_groups_table)
combine_perturbations <- function(perturbation_matrix, gRNA_groups_table) {
  # check that gRNA_id and gRNA_group are present in `gRNA_groups` df
  if (!all(c("gRNA_id", "gRNA_group") %in% colnames(gRNA_groups_table))) {
    stop("`gRNA_id` and `gRNA_group` should be columns of `gRNA_groups`.")
  }

  # convert the data frame to a named list
  gRNA_groups <- as.character(unique(gRNA_groups_table$gRNA_group))
  gRNA_groups <- magrittr::set_names(gRNA_groups, gRNA_groups)
  gRNA_grp_list <- lapply(X = gRNA_groups, FUN = function(gRNA_group) {
    dplyr::filter(gRNA_groups_table, gRNA_group == !!gRNA_group) %>% dplyr::pull(gRNA_id)
  })
  # Determine which gRNAs are not contained within gRNA_grps
  leftover_gRNAs <- setdiff(row.names(perturbation_matrix), unlist(gRNA_grp_list))
  out_leftover <- perturbation_matrix[leftover_gRNAs,]
  out_grped <- sapply(X = names(gRNA_grp_list), FUN = function(grp_name) {
    mat_sub <- perturbation_matrix[gRNA_grp_list[[grp_name]],]
    Matrix::colSums(mat_sub)
  }) %>% Matrix::t()
  out_grped <- out_grped >= 1
  out <- rbind(out_leftover, out_grped)
}


#' Threshold a gRNA count matrix
#'
#' Thresholds a given gRNA count matrix.
#'
#' See "Exponential family measurement error models for single-cell CRISPR screens" by Barry et al. 2022 for an in-depth analysis of thresholding in single-cell CRISPR screens.
#'
#' @param gRNA_matrix a gRNA-by-cell expression matrix; the matrix can be represented as a sparse matrix (as implemented by the Matrix package) or a dense matrix (as implemented by base R)
#' @param threshold the threshold used to assign perturbations to cells; counts above the threshold are set to 1 (indicating "perturbed"), and counts below the threshold are set to 0 (indicating "unperturbed")
#'
#' @return a binary matrix of imputed perturbation assignments
#' @export
#'
#' @examples
#' data(gRNA_matrix)
#' perturbation_matrix <- threshold_gRNA_matrix(gRNA_matrix)
threshold_gRNA_matrix <- function(gRNA_matrix, threshold = 3) {
  return(gRNA_matrix >= threshold)
}
