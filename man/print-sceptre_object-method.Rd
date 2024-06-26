% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s4_helpers.R
\name{print,sceptre_object-method}
\alias{print,sceptre_object-method}
\title{Print}
\usage{
\S4method{print}{sceptre_object}(x)
}
\arguments{
\item{x}{a \code{sceptre_object}}
}
\value{
the value NULL
}
\description{
\code{print()} prints information about the dataset and the status of the analysis to the console. The output contains several fields: Attributes of the data (summarizing key features of the data), Analysis status (indicating the analysis functions that have been called), Analysis parameters (summarizing the analysis parameters set in \code{set_analysis_parameters()}), gRNA-to-cell assignment information (summarizing the outcome of the gRNA-to-cell assignment step), and Summary of results (summarizing the key analysis results). A subset of these fields may be printed, depending on the status of the analysis.
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
print(sceptre_object)
}
