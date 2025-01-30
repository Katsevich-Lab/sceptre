#' Gene position data frames
#'
#' `gene_position_data_frame_grch38` and `gene_position_data_frame_grch37` contain the coordinate and transcription start site position of each gene relative to reference genome GRCh38 and GRCh37, respectively. Both `gene_position_data_frame_grch38` and `gene_position_data_frame_grch37` were constructed from reference genomes available on the 10x Genomics website. The GRCh38 reference genome has been used by 10x Cell Ranger since 2020.
#'
#' @usage data(gene_position_data_frame_grch38)
#' @examples
#' head(gene_position_data_frame_grch38)
#' head(gene_position_data_frame_grch37)
"gene_position_data_frame_grch38"

#' @rdname gene_position_data_frame_grch38
#' @usage data(gene_position_data_frame_grch37)
"gene_position_data_frame_grch37"

#' Example high-MOI data
#'
#' The example high-MOI CRISPRi data are a small simulated dataset modeled on
#' that of "A genome-wide framework for mapping gene regulation via cellular
#' genetic screens" by Gasperini et al., 2019. These data contain perturbations
#' both of gene transcription start sites (TSSs) and enhancers, as well as non-
#' targeting perturbations.
#'
#' The example data are stored in a list containing the following components:
#' - `response_matrix`: the gene-by-cell expression matrix
#' - `grna_matrix`: the gRNA-by-cell expression matrix
#' - `extra_covariates`: a data frame containing a single column, `batch`, specifying the batch in which each cell was sequenced (`batch_1` or `batch_2`)
#' - `gene_names`: the human-readable name of each gene
#'
#' @usage data(highmoi_example_data)
"highmoi_example_data"

#' gRNA target data frame
#'
#' The gRNA target data frame corresponding to the example high MOI data. The data frame contains the columns `grna_id` (ID of an individual gRNA), `grna_target` (genomic target of the gRNA), `chr` (target chromosome), `start` (start coordinate of target location), and `end` (end coordinate of target location).
#' @usage data(grna_target_data_frame_highmoi)
"grna_target_data_frame_highmoi"

#' Example low-MOI data
#'
#' The example low-MOI data are a small simulated dataset modeled on that in the
#' paper "Characterizing the molecular regulation of inhibitory immune
#' checkpoints with multimodal single-cell screens" by Papalexi et al., 2021.
#' This dataset includes gene-targeting CRISPRko perturbations as well as non-
#' targeting perturbations.
#'
#' The example data are stored in a list containing four components:
#' - `response_matrix`: the gene-by-cell expression matrix
#' - `grna_matrix`: the gRNA-by-cell expression matrix
#' - `grna_target_data_frame`:  a data frame containing the columns `grna_id` (ID of an individual gRNA) and `grna_target` (genomic target of the gRNA)
#' - `extra_covariates`: a data frame with a single column, `batch`, specifying the batch in which each cell was sequenced (`batch_1` or `batch_2`)
#'
#' @usage data(lowmoi_example_data)
"lowmoi_example_data"
