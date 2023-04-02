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
#' @usage data(gene_matrix_highmoi)
"gene_matrix_highmoi"

#' gRNA expression matrix
#'
#' A gRNA-by-cell expression matrix.
#'
#' @inheritSection highmoi_data Overview
#' @usage data(gRNA_matrix_highmoi)
"gRNA_matrix_highmoi"

#' Covariate matrix
#'
#' A matrix of cell-specific technical factors:
#' \describe{
#'   \item{`lg_gRNA_lib_size`}{(numeric) log-transformed gRNA library size}
#'   \item{`lg_gene_lib_size`}{(numeric) log-transformed gene library size}
#'   \item{`p_mito`}{(numeric) fraction of sequenced gene transcripts that map to mitochondrial genes}
#'   \item{`batch`}{(factor) sequencing batch (`batch_1`, `batch_2`)}
#' }
#' @inheritSection highmoi_data Overview
#' @usage data(covariate_matrix_highmoi)
"covariate_matrix_highmoi"

#' Gene-gRNA group pairs
#'
#' The pairs of genes and gRNA groups that we seek to test for association. Columns `gene_id` (required), `gRNA_group` (required), and `pair_type` (optional).
#' @inheritSection highmoi_data Overview
#' @usage data(gene_gRNA_group_pairs_highmoi)
"gene_gRNA_group_pairs_highmoi"

#' gRNA groups table
#'
#' A data frame that maps each gRNA ID to its gRNA group and gRNA type. The columns are `gRNA_id` (required), `gRNA_group` (required), and `gRNA_type` (optional).
#' @inheritSection highmoi_data Overview
#' @usage data(gRNA_groups_table_highmoi)
"gRNA_groups_table_highmoi"
