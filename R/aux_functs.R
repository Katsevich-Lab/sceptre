#' Combine gRNAs (deprecated)
#'
#' Combines gRNAs that target the same site into a single "combined" gRNA via addition of the expression levels.
#'
#' Combining gRNAs that target the same site helps to increase statistical power. This operation is especially appealing when gRNA expression levels are low or there are not a lot of (e.g., < 100,000) cells.
#'
#' @param gRNA_matrix a gRNA-by-cell expression matrix; the row names should be gRNA IDs
#' @param site_table a data frame with column names "site" and "gRNA_id". The "site" column should give the target site of each gRNA ID
#' @return a "collapsed" gRNA-by-cell expression matrix, where gRNAs that target the same site have been collapsed into a single row (via addition of the constituent gRNAs)
#'
#' @examples
#' \dontrun{
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
#' }
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


#' Create perturbation matrix
#'
#' A threshold of 3 is a reasonable default for high-MOI single-cell CRISPR screens; see "Exponential family measurement error models for single-cell CRISPR screens" by Barry et al. 2022 for an in-depth analysis of thresholding and more sophisticated thresholding techniques.
#'
#' @param gRNA_matrix a gRNA-by-cell expression matrix stored as a sparse matrix (as implemented by the Matrix package) or a dense matrix (as implemented by base R); the row names should be the gRNA IDs
#' @param gRNA_groups_table a data frame with columns `gRNA_id` and `gRNA_group`.
#' @param threshold
#'
#' @return
#' @export
#'
#' @examples
#' data("gRNA_matrix")
#' data("gRNA_groups")
#' gRNA_groups_table <- gRNA_groups
create_perturbation_matrix <- function(gRNA_matrix, gRNA_groups_table, threshold = 3) {
  # check that gRNA_id and gRNA_group are present in `gRNA_groups` df
  if (!all(c("gRNA_id", "gRNA_group") %in% colnames(gRNA_groups_table))) {
    stop("`gRNA_id` and `gRNA_group` should be columns of `gRNA_groups`.")
  }
  # first, threshold the gRNAs, thereby assigning perturbation indicators
  gRNA_matrix_thresh <- gRNA_matrix >= threshold

  # convert the data frame to a named list
  gRNA_groups <- as.character(unique(gRNA_groups_table$gRNA_group))
  gRNA_groups <- magrittr::set_names(gRNA_groups, gRNA_groups)
  gRNA_grp_list <- lapply(X = gRNA_groups, FUN = function(gRNA_group) {
    dplyr::filter(gRNA_groups_table, gRNA_group == !!gRNA_group) %>% dplyr::pull(gRNA_id)
  })
  # Determine which gRNAs are not contained within gRNA_grps
  leftover_gRNAs <- setdiff(row.names(gRNA_matrix), unlist(gRNA_grp_list))
  out_leftover <- gRNA_matrix[leftover_gRNAs,]
  out_grped <- sapply(X = names(gRNA_grp_list), FUN = function(grp_name) {
    mat_sub <- gRNA_matrix_thresh[gRNA_grp_list[[grp_name]],]
    Matrix::colSums(mat_sub)
  }) %>% Matrix::t()
  out_grped <- out_grped >= 1
  out <- rbind(out_leftover, out_grped)
}


#' Threshold a gRNA count matrix
#'
#' Thresholds a given gRNA count matrix
#'
#' @param gRNA_matrix a gRNA-by-cell expression matrix; the matrix can be represented as a sparse matrix (as implemented by the Matrix package) or a dense matrix (as implemented by base R)
#' @param threshold the threshold used to assign perturbations to cells; counts above the threshold are set to 1 (indicating "perturbed"), and counts below the threshold are set to 0 (indicating "unperturbed")
#'
#' @return a binary matrix of imputed perturbation assignments
#' @export
#'
#' @examples
#' data(gRNA_matrix)
#' gRNA_matrix_binary <- threshold_gRNA_matrix(gRNA_matrix)
threshold_gRNA_matrix <- function(gRNA_matrix, threshold = 3) {
  return(gRNA_matrix >= threshold)
}

