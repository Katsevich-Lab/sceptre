% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plot_run_calibration_check}
\alias{plot_run_calibration_check}
\title{Plot run calibration check}
\usage{
plot_run_calibration_check(
  sceptre_object,
  point_size = 0.55,
  transparency = 0.8,
  return_indiv_plots = FALSE
)
}
\arguments{
\item{sceptre_object}{a \code{sceptre_object} that has had \code{run_calibration_check} called on it}

\item{point_size}{(optional; default \code{0.55}) the size of the individual points in the plot}

\item{transparency}{(optional; default \code{0.8}) the transparency of the individual points in the plot}

\item{return_indiv_plots}{(optional; default \code{FALSE}) if \code{FALSE} then a list of \code{ggplot} is returned; if \code{TRUE} then a single \code{cowplot} object is returned.}
}
\value{
a single \code{cowplot} object containing the combined panels (if \code{return_indiv_plots} is set to \code{TRUE}) or a list of the individual panels (if \code{return_indiv_plots} is set to \code{FALSE})
}
\description{
\code{plot_run_calibration_check()} creates a visualization of the outcome of the calibration check. The visualization consists of four panels, which we describe below.
\itemize{
\item The upper left panel is a QQ plot of the p-values plotted on an untransformed scale. The p-values ideally should lie along the diagonal line, indicating uniformity of the p-values in the \emph{bulk} of the distribution.
\item The upper right panel is a QQ plot of the p-values plotted on a negative log-10 transformed scale. The p-values ideally should lie along the diagonal line (with the majority of the p-values falling within the gray confidence band), indicating uniformity of the p-values in the \emph{tail} of the distribution.
\item The lower left panel is a histogram of the estimated log-2 fold changes. The histogram ideally should be roughly symmetric and centered around zero.
\item Finally, the bottom right panel is a text box displaying (i) the number of false discoveries that \code{sceptre} has made on the negative control data and (ii) the mean estimated log-fold change.
}
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
sceptre_object |>
  set_analysis_parameters(
    side = "left",
    resampling_mechanism = "permutations"
  ) |>
  assign_grnas(method = "thresholding") |>
  run_qc() |>
  run_calibration_check(
    parallel = TRUE,
    n_processors = 2,
    n_calibration_pairs = 500,
    calibration_group_size = 2,
  ) |>
  plot_run_calibration_check()
}
