% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_functs.R
\name{import_data_from_parse}
\alias{import_data_from_parse}
\title{Import data from Parse (experimental)}
\usage{
import_data_from_parse(
  gene_mat_fp,
  grna_mat_fp,
  all_genes_fp,
  all_grnas_fp,
  moi,
  grna_target_data_frame,
  extra_covariates = data.frame()
)
}
\arguments{
\item{gene_mat_fp}{file path to the gene \code{count_matrix.mtx} file.}

\item{grna_mat_fp}{file path to the gRNA \code{count_matrix.mtx} file.}

\item{all_genes_fp}{file path to the \code{all_genes.csv} file containing the gene IDs.}

\item{all_grnas_fp}{file path to the \code{all_guides.csv} file containing the gRNA IDs. The gRNA IDs are assumed to be in the second column (i.e., the "gene_name" column) of this file.}

\item{moi}{a string indicating the MOI of the dataset, either "low" or "high".}

\item{grna_target_data_frame}{a data frame containing columns \code{grna_id} and \code{grna_target} mapping each individual gRNA to its target. Non-targeting gRNAs should be assigned a label of "non-targeting". Optionally, \code{grna_target_data_frame} can contain columns \code{chr}, \code{start}, and \code{end}, giving the chromosome, start coordinate, and end coordiante, respectively, of each gRNA. Additionally, \code{grna_target_data_frame} can contain the column \code{vector_id} specifying the vector to which a given gRNA belongs.}

\item{extra_covariates}{(optional) a data frame containing extra covariates (e.g., batch, biological replicate) beyond those that \code{sceptre} can compute.}
}
\value{
an initialized \code{sceptre_object}
}
\description{
\code{import_data_from_parse()} imports data from the output of the Parse count matrix generation program. See \href{https://timothy-barry.github.io/sceptre-book/import-data.html#import-from-the-parse-program-experimental}{Chapter 1 of the manual} for more information about this function.
}
\details{
\code{import_data_from_parse()} is experimental, and the API of this function is subject to change. We expect the API to solidify as we learn more about the Parse platform and the structure of the Parse count matrix generation program output.
}
\examples{
directory <- paste0(
  system.file("extdata", package = "sceptredata"),
  "/parse_example/"
)
gene_mat_fp <- paste0(directory, "gene_mat.mtx")
grna_mat_fp <- paste0(directory, "grna_mat.mtx")
all_genes_fp <- paste0(directory, "all_genes.csv")
all_grnas_fp <- paste0(directory, "all_grnas.csv")
grna_target_data_frame <- data.frame(
  grna_id = c("guide_A", "guide_B", "guide_C"),
  grna_target = c("target-A", "target-B", "non-targeting")
)
sceptre_object <- import_data_from_parse(
  gene_mat_fp = gene_mat_fp,
  grna_mat_fp = grna_mat_fp,
  all_genes_fp = all_genes_fp,
  all_grnas_fp = all_grnas_fp,
  moi = "low",
  grna_target_data_frame = grna_target_data_frame
)
}
