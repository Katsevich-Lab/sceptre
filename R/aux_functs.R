#' Combine gRNAs
#'
#' Combines gRNAs that target the same site into a single "combined" gRNA via addition of the constituent expression levels.
#'
#' It typically is a good idea to combine gRNAs that target the same site into a single site-specific "combined" gRNA to increase statistical power. This is especially the case when gRNA expression levels are low or there are not a lot of (e.g., < 100,000) cells.
#'
#' @param gRNA_matrix a gRNA-by-cell expression matrix; the matrix should have row names giving the gRNA IDs
#' @param gRNA_grps a named list; each entry of the list should be a character vector giving the IDs of gRNAs that target the same site, with the name of that entry being the name of the targeted site
#' @return a "combined" gRNA-by-cell expression matrix, where gRNAs that target the same site have been collapsed into a single row (via addition of the constituent gRNAs)
#' @export
#'
#' @examples
#' library(magrittr)
#' data(gRNA_matrix)
#' data(gRNA_grps)
#' combined_gRNA_matrix <- combine_gRNAs(gRNA_matrix, gRNA_grps)
#'
#' # Example explanation:
#' # `gRNA_grps` indicates that the gRNAs "AGAAGAGTCAGCCTTGCAC,"
#' # "GAAGAGTCAGCCTTGCACT," "GTTTAGGGAACCCAGTGCA," "GTAACTTCATTTGCAGCAA,"
#' # "TTACTTTTTATCAAGCCAA," "TCTAATTTAAGACCTGGGT," "AAGGTATCATTTTCTTGTC,"
#' # "TTAGCCATAGAACCACTAC," "TTTTCCAGTAGTGGTTCTA," and "TACATGGCATAGATCTTTA"
#' # all target the site "chr8:128428069-128428469". The function  collapses
#' # the above gRNAs into a single "combined" gRNA called
#' # "chr8:128428069-128428469" (and similarly for the other sites).
combine_gRNAs <- function(gRNA_matrix, gRNA_grps) {
  out <- sapply(X = names(gRNA_grps), FUN = function(grp_name) {
    # cat(paste0("Operating on ", grp_name, ".\n"))
    mat_sub <- gRNA_matrix[gRNA_grps[[grp_name]],]
    Matrix::colSums(mat_sub)
  }) %>% Matrix::t()
  return(out)
}
