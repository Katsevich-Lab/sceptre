% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pair_constructor_functs.R
\name{construct_cis_pairs}
\alias{construct_cis_pairs}
\title{Construct \emph{cis} pairs}
\usage{
construct_cis_pairs(
  sceptre_object,
  positive_control_pairs = data.frame(),
  distance_threshold = 500000L,
  response_position_data_frame = gene_position_data_frame_grch38
)
}
\arguments{
\item{sceptre_object}{a \code{sceptre_object}}

\item{positive_control_pairs}{(optional) a data frame with columns \code{grna_target} and \code{response_id} containing the positive control pairs; if supplied, the positive control targets are excluded from the \emph{cis} pairs.}

\item{distance_threshold}{(optional) target-response pairs located within \code{distance_threshold} bases of one another and on the same chromosome are included in the \emph{cis} discovery set.}

\item{response_position_data_frame}{(optional) a data frame with columns \code{response_id}, \code{chr}, and \code{position} giving the genomic coordinate of each response; by default \code{response_position_data_frame} is set to a data frame containing the genomic coordinate of each gene in the human genome relative to reference genome GRCh38.}
}
\value{
a data frame with columns \code{grna_target} and \code{response_id} containing the \emph{cis} discovery set
}
\description{
\code{construct_cis_pairs()} is a helper function to facilitate construction the \emph{cis} pairs. See \href{https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html#sec-set-analysis-parameters_construct_cis_pairs}{Section 2.2.2} of the manual for more detailed information about this function.
}
