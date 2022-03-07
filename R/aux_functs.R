#' Combine gRNAs
#'
#'
#'
#' @param gRNA_matrix XX
#' @param gRNA_grps XX
#'
#' @return
#' @export
#'
#' @examples
#' library(magrittr)
#' data(gRNA_matrix)
#' indiv_gRNAs <- sort(row.names(gRNA_matrix))
#' grp_size <- 5L
#' idxs <- seq(from = 1, to = length(indiv_gRNAs) + grp_size, by = grp_size)
#' group_idxs <- seq(1, length(idxs) - 1)
#' gRNA_grps <- lapply(X = group_idxs, FUN = function(i) {
#'  gRNA_grp <- indiv_gRNAs[seq(idxs[i], idxs[i + 1] - 1)]
#' }) %>% set_names(., paste0("gRNA_grp_",  group_idxs))
#' combined_gRNA_matrix <- combine_gRNAs(gRNA_matrix, gRNA_grps)
combine_gRNAs <- function(gRNA_matrix, gRNA_grps) {
  out <- sapply(X = names(gRNA_grps), FUN = function(grp_name) {
    # cat(paste0("Operating on ", grp_name, ".\n"))
    mat_sub <- gRNA_matrix[gRNA_grps[[grp_name]],]
    Matrix::colSums(mat_sub)
  }) %>% Matrix::t()
  return(out)
}
