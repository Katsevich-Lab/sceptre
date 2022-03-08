#' `sceptre` package example data
#' @name data
#' @section Overview:
#' This object is part of the example data in the `sceptre` package. The data are taken from the paper "Global Analysis of Enhancer Targets Reveals Convergent Enhancer-Driven Regulatory Modules" by Xie et al., 2019. The authors used CRISPRi-based perturb-seq to target >500 putative enhancers in a population of K562 cells. The example data consist of a subset of the sequenced gRNAs and genes measured on 106,660 cells.
NULL


#' Gene expression matrix
#'
#' A gene-by-cell expression matrix.
#'
#' @inheritSection data Overview
#' @usage data(gene_matrix)
"gene_matrix"

#' gRNA expression matrix
#'
#' A gRNA-by-cell expression matrix.
#'
#' @inheritSection data Overview
#' @usage data(gRNA_matrix)
"gRNA_matrix"

#' Covariate matrix
#'
#' A matrix of cell-specific technical factors:
#' \describe{
#'   \item{`lg_gene_lib_size`}{(numeric) log-transformed gene library size}
#'   \item{`lg_gRNA_lib_size`}{(numeric) log-transformed gRNA library size}
#'   \item{`batch`}{(factor) sequencing batch (batch_1, ..., batch_5)}
#'   \item{`p_mito`}{(numeric) fraction of sequenced gene transcripts that map to mitochondrial genes}
#' }
#' @inheritSection data Overview
#' @usage data(covariate_matrix)
"covariate_matrix"

#' Gene-gRNA pairs
#'
#' The gene-gRNA pairs that we seek to test for association. Columns include `gene_id` and `gRNA_id`.
#' @inheritSection data Overview
#' @usage data(gene_gRNA_pairs)
"gene_gRNA_pairs"

#' gRNA groups
#'
#' A named list grouping gRNAs by the chromosomal site that they target.
#' @inheritSection data Overview
#' @usage data(gRNA_grps)
"gRNA_grps"
