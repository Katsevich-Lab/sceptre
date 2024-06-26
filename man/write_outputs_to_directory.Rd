% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s4_analysis_functs_2.R
\name{write_outputs_to_directory}
\alias{write_outputs_to_directory}
\title{Write outputs to directory}
\usage{
write_outputs_to_directory(sceptre_object, directory)
}
\arguments{
\item{sceptre_object}{a \code{sceptre_object}}

\item{directory}{a string giving the file path to a directory on disk in which to write the results}
}
\value{
the value NULL
}
\description{
\code{write_outputs_to_directory()} writes the outputs of a \code{sceptre} analysis to a directory on disk. \code{write_outputs_to_directory()} writes several files to the specified directory: a text-based summary of the analysis (\code{analysis_summary.txt}), the various plots (\verb{*.png}), the calibration check, power check, discovery analysis results (\code{results_run_calibration_check.rds}, \code{results_run_power_check.rds}, and \code{results_run_discovery_analysis.rds}, respectively), and the binary gRNA-to-cell assignment matrix (\code{grna_assignment_matrix.rds}). See \href{https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_write_outputs_to_directory}{Section 8 of the introductory chapter in the manual} for more information about this function.
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
positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
discovery_pairs <- construct_cis_pairs(sceptre_object,
  positive_control_pairs = positive_control_pairs,
  distance_threshold = 5e6
)
sceptre_object |>
  set_analysis_parameters(
    side = "left",
    resampling_mechanism = "permutations",
    discovery_pairs = discovery_pairs,
    positive_control_pairs = positive_control_pairs
  ) |>
  assign_grnas(method = "thresholding") |>
  run_qc() |>
  run_calibration_check(
    parallel = TRUE,
    n_processors = 2
  ) |>
  run_power_check() |>
  run_discovery_analysis(
    parallel = TRUE,
    n_processors = 2
  ) |>
  write_outputs_to_directory(paste0(tempdir(), "/sceptre_outputs"))
# files written to "sceptre_outputs" in tempdir()
}
