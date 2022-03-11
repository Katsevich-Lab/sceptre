#' Combine gRNAs
#'
#' Combines gRNAs that target the same site into a single "combined" gRNA via addition of the expression levels.
#'
#' Combining gRNAs that target the same site helps to increase statistical power. This operation is especially appealing when gRNA expression levels are low or there are not a lot of (e.g., < 100,000) cells.
#'
#' @param gRNA_matrix a gRNA-by-cell expression matrix; the row names should be gRNA IDs
#' @param site_table a data frame with column names "site" and "gRNA_id". The "site" column should give the target site of each gRNA ID
#' @return a "collapsed" gRNA-by-cell expression matrix, where gRNAs that target the same site have been collapsed into a single row (via addition of the constituent gRNAs)
#' @export
#'
#' @examples
#' library(magrittr)
#' library(Matrix)
#' library(dplyr)
#'
#' # 1. First example
#' data(gRNA_matrix)
#' data(site_table)
#' combined_gRNA_matrix <- combine_gRNAs(gRNA_matrix, site_table)
#'
#' # 2. Second example
#' # Here we group only the gRNAs that target sites "chr10:17457016-17457416"
#' # and "chr18:48566684-48567084;" all others remain ungrouped
#' site_table_2 <- site_table %>%
#' filter(site %in% c("chr10:17457016-17457416", "chr18:48566684-48567084"))
#' combined_gRNA_matrix_2 <- combine_gRNAs(gRNA_matrix, site_table_2)
combine_gRNAs <- function(gRNA_matrix, site_table) {
  # check that "site" and "gRNA_id" are present as column names
  if (!all(c("site", "gRNA_id") %in% colnames(site_table))) {
    stop("`site` and `gRNA_id` should be columns of `site_table`.")
  }
  # convert the data frame to a named list
  sites <- unique(site_table$site)
  sites <- magrittr::set_names(sites, sites)
  gRNA_grps <- lapply(X = sites, FUN = function(site) {
    dplyr::filter(site_table, site == !!site) %>% dplyr::pull(gRNA_id)
  })
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
