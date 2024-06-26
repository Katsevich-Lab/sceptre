% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plot_response_grna_target_pair}
\alias{plot_response_grna_target_pair}
\title{Plot response-gRNA-target pair}
\usage{
plot_response_grna_target_pair(sceptre_object, response_id, grna_target)
}
\arguments{
\item{sceptre_object}{a \code{sceptre_object} that has had \code{run_qc()} called on it}

\item{response_id}{a string containing a response ID}

\item{grna_target}{a string containing a gRNA target (or, if \code{grna_integration_strategy} is set to \code{"singleton"}, an individual gRNA ID)}
}
\value{
a violin plot
}
\description{
\code{plot_response_grna_target_pair()} creates a violin plot of the expression level of a given response as a function of the "treatment status" (i.e., treatment or control) of a given gRNA target. The left (resp., right) violin plot shows the expression level of the response in treatment (resp., control) cells. The expression level is normalized by dividing by \code{n_response_umis}, adding a pseudo-count of 1 and then taking the log transform. If the given response-gRNA-target pair has been analyzed, the p-value for the test of association also is displayed.
}
\details{
If \code{grna_integration_strategy} is set to \code{"singleton"}, then \code{grna_target} should be set to a gRNA ID.
}
\examples{
library(sceptredata)
data(highmoi_example_data)
data(grna_target_data_frame_highmoi)
# import data
sceptre_object <- import_data(
  response_matrix = highmoi_example_data$response_matrix,
  grna_matrix = highmoi_example_data$grna_matrix,
  grna_target_data_frame = grna_target_data_frame_highmoi,
  moi = "high",
  extra_covariates = highmoi_example_data$extra_covariates,
  response_names = highmoi_example_data$gene_names
)
discovery_pairs <- construct_cis_pairs(sceptre_object)
sceptre_object |>
  set_analysis_parameters(
    side = "left",
    discovery_pairs = discovery_pairs,
    resampling_mechanism = "permutations",
  ) |>
  assign_grnas(method = "thresholding") |>
  run_qc() |>
  run_discovery_analysis(
    parallel = TRUE,
    n_processors = 2
  ) |>
  plot_response_grna_target_pair(
    response_id = "ENSG00000136938",
    grna_target = "candidate_enh_20"
  )
}
