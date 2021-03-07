#' Single-cell CRISPR screen example data
#'
#' Example single-cell CRISPR screen data from the study A Genome-wide Framework for Mapping Gene Regulation via Cellular Genetic Screens by Gasperini et al. The authors used CRISPRi to perturb 6,000 candidate enhancers and assessed the impact of these perturbations via scRNA-seq. The example data come from a single gRNA (ID: chr7.3255_top_two), a single gene (ID: ENSG00000164713), and the set of cell-specific technical factors. There are 205,797 cells in total.
#' @name data
#' @format
#' - expressions: an integer vector containing the gene expression data across all cells, measured in UMIs.

#' - gRNA_indicators: a binary (i.e., 0/1) vector indicating whether the perturbation was detected in a given cell.
#'
#' - covariate_matrix: the set of cell-specific technical factors. The following columns are in the data frame:
#' \describe{
#'   \item{p_mito}{fraction of UMIs that mapped to mitochondrial genes}
#'   \item{prep_batch}{sequencing batch (either prep_batch_1 or prep_batch_2)}
#'   \item{lg_total_umis}{log-transformed total UMI count}
#'   \item{lg_guide_count}{log-transformed number of perturbations detected}
#'   \item{lg_n_genes}{log-transformed number of genes expressed}
#' }
#' @source A genome-wide framework for Mapping Gene Regulation via Cellular Genetic Screens, Gasperini et al. 2019. Link: \url{https://www.sciencedirect.com/science/article/pii/S009286741831554X?via%3Dihub}
NULL

#' @rdname data
"expressions"

#' @rdname data
"gRNA_indicators"

#' @rdname data
"covariate_matrix"
