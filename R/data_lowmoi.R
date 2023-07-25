#' \code{sceptre} package example low MOI data
#'
#' @name lowmoi_data
#' @section Overview:
#' The example low MOI data in the \code{sceptre} package are taken from the paper "Characterizing the molecular regulation of inhibitory immune checkpoints with multimodal single-cell screens" by Papalexi et al., 2021. The authors used CRISPRko to target 26 genes in 20,729 THP-1 cells and measured the effects of these perturbations via single-cell sequencing. The authors designed three or four gRNAs to target each gene, as well as a library of nine non-targeting gRNAs. The example data include 290 (or 2\%) randomly selected genes from the original dataset.
NULL


#' Response expression matrix
#'
#' A gene-by-cell expression matrix. The row names of the matrix are the names of the genes.
#'
#' @inheritSection lowmoi_data Overview
#' @usage data(response_matrix_lowmoi)
"response_matrix_lowmoi"


#' gRNA expression matrix
#'
#' A gRNA-by-cell expression matrix. The row names of the matrix are the names of the individual gRNAs.
#'
#' @inheritSection lowmoi_data Overview
#' @usage data(grna_matrix_lowmoi)
"grna_matrix_lowmoi"


#' Extra covariates
#'
#' A data frame containing the cell-specific technical factors beyond those that `sceptre` can compute:
#' \describe{
#'   \item{\code{bio_rep}}{(factor) the biological replicate in which the cell was sequenced (\code{rep_1}, \code{rep_2}, or \code{rep_3})}
#' }
#' @inheritSection lowmoi_data Overview
#' @usage data(covariate_data_frame_lowmoi)
"extra_covariates_lowmoi"


#' gRNA group data frame
#'
#' A data frame containing the gRNA group information. The data frame has two columns: \code{grna_id} and \code{grna_group}. The former column provides the ID of each individual gRNA, and the latter specifies the group to which each individual gRNA belongs. gRNAs that target the same gene transcription start site are assigned to the same group. Non-targeting gRNAs are assigned to a group labeled "non-targeting".
#'
#' @inheritSection lowmoi_data Overview
#' @usage data(grna_group_data_frame_lowmoi)
"grna_group_data_frame_lowmoi"
