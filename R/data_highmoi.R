#' \code{sceptre} package example high MOI data
#' @name highmoi_data
#' @section Overview:
#' This object is part of the example data. The data are taken from the paper "A genome-wide framework for mapping gene regulation via cellular genetic screens" by Gasperini et al., 2019. The authors used a CRISPRi-based assay to target 5,000+ putative enhancers in a population of K562 cells. The authors additionally targeted 200+ gene transcription start sites (TSSs) and designed a library of 50 non-targeting gRNAs to serve as negative controls. Genes, gRNAs, and cells are all down-sampled to reduce the size of the data.
NULL


#' Gene expression matrix
#'
#' A gene-by-cell expression matrix.
#'
#' @inheritSection highmoi_data Overview
#' @usage data(response_matrix_highmoi_experimental)
"response_matrix_highmoi_experimental"

#' gRNA expression matrix
#'
#' A gRNA-by-cell expression matrix.
#'
#' @inheritSection highmoi_data Overview
#' @usage data(grna_matrix_highmoi_experimental)
"grna_matrix_highmoi_experimental"

#' Covariate matrix
#'
#' A matrix of cell-specific technical factors:
#' \describe{
#'   \item{`batch`}{(factor) sequencing batch (`batch_1`, `batch_2`)}
#' }
#' @inheritSection highmoi_data Overview
#' @usage data(extra_covariates_highmoi_experimental)
"extra_covariates_highmoi_experimental"

#' Discovery gene-gRNA group pairs
#'
#' The discovery response-gRNA group pairs that we seek to test for association, with columns `response_id` (required), `grna_group` (required), and `type` (optional).
#' @inheritSection highmoi_data Overview
#' @usage data(discovery_pairs_highmoi_experimental)
"discovery_pairs_highmoi_experimental"

#' Positive control gene-gRNA group pairs
#'
#' The positive control response-gRNA group pairs, with columns `response_id` (required), `grna_group` (required), and `type` (optional).
#' @inheritSection highmoi_data Overview
#' @usage data(discovery_pairs_highmoi_experimental)
"pc_pairs_highmoi_experimental"

#' gRNA groups table
#'
#' A data frame that maps each gRNA ID to its gRNA group.
#' @inheritSection highmoi_data Overview
#' @usage data(grna_group_data_frame_highmoi_experimental)
"grna_group_data_frame_highmoi_experimental"
