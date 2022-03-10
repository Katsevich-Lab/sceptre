#' Combine gRNAs
#'
#' Combines gRNAs that target the same site into a single "combined" gRNA via addition of the expression levels.
#'
#' It typically is a good idea to combine gRNAs that target the same site into a single site-specific "combined" gRNA to increase statistical power. This is especially true when gRNA expression levels are low or there are not a lot of (e.g., < 100,000) cells.
#'
#' @param gRNA_matrix a gRNA-by-cell expression matrix; the row names should be gRNA IDs
#' @param gRNA_grps a named list; each entry of the list should be a character vector giving the IDs of gRNAs that target the same site, with the name of that entry being the name of the targeted site
#' @return a "combined" gRNA-by-cell expression matrix, where gRNAs that target the same site have been collapsed into a single row (via addition of the constituent gRNAs)
#' @export
#'
#' @examples
#'
#' # 1. First example
#' # `gRNA_grps` indicates that the gRNAs "AGAAGAGTCAGCCTTGCAC,"
#' # "GAAGAGTCAGCCTTGCACT," "GTTTAGGGAACCCAGTGCA," "GTAACTTCATTTGCAGCAA,"
#' # "TTACTTTTTATCAAGCCAA," "TCTAATTTAAGACCTGGGT," "AAGGTATCATTTTCTTGTC,"
#' # "TTAGCCATAGAACCACTAC," "TTTTCCAGTAGTGGTTCTA," and "TACATGGCATAGATCTTTA"
#' # all target the site "chr8:128428069-128428469". The function  collapses
#' # the above gRNAs into a single "combined" gRNA called
#' # "chr8:128428069-128428469" (and similarly for the other sites).
#' library(magrittr)
#' library(Matrix)
#'
#' data(gRNA_matrix)
#' data(gRNA_grps)
#' combined_gRNA_matrix <- combine_gRNAs(gRNA_matrix, gRNA_grps)
#'
#' # 2. Second example
#' # Here we group only the gRNAs in the first two gRNA groups;
#' # all other gRNAs remain ungrouped.
#' gRNA_grps <- gRNA_grps[1:2]
#' combined_gRNA_matrix_2 <- combine_gRNAs(gRNA_matrix, gRNA_grps)
combine_gRNAs <- function(gRNA_matrix, gRNA_grps) {
  # Determine which gRNAs are not contained within gRNA_grps
  leftover_gRNAs <- setdiff(row.names(gRNA_matrix), unlist(gRNA_grps))
  out_leftover <- gRNA_matrix[leftover_gRNAs,]
  out_grped <- sapply(X = names(gRNA_grps), FUN = function(grp_name) {
    mat_sub <- gRNA_matrix[gRNA_grps[[grp_name]],]
    Matrix::colSums(mat_sub)
  }) %>% Matrix::t()
  out <- rbind(out_leftover, out_grped)
  return(out)
}

