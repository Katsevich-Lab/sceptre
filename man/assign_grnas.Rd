% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s4_analysis_functs_1.R
\name{assign_grnas}
\alias{assign_grnas}
\title{Assign gRNAs to cells}
\usage{
assign_grnas(
  sceptre_object,
  method = "default",
  print_progress = TRUE,
  parallel = FALSE,
  n_processors = "auto",
  log_dir = tempdir(),
  ...
)
}
\arguments{
\item{sceptre_object}{a \code{sceptre_object}}

\item{method}{(optional) a string indicating the method to use to assign the gRNAs to cells, one of \code{"mixture"}, \code{"thresholding"}, or \code{"maximum"}}

\item{print_progress}{(optional; default \code{TRUE}) a logical indicating whether to print progress updates}

\item{parallel}{(optional; default \code{FALSE}) a logical indicating whether to run the function in parallel}

\item{n_processors}{(optional; default "auto") an integer specifying the number of processors to use if \code{parallel} is set to \code{TRUE}. The default, \code{"auto"}, automatically detects the number of processors available on the machine.}

\item{log_dir}{(optional; default \code{tempdir()}) a string indicating the directory in which to write the log files (ignored if \code{parallel = FALSE})}

\item{...}{optional method-specific additional arguments}
}
\value{
an updated \code{sceptre_object} in which the gRNA assignments have been carried out
}
\description{
\code{assign_grnas()} performs the gRNA-to-cell assignments. \code{sceptre} provides three gRNA-to-cell assignment strategies: the mixture method, the thresholding method, and the maximum method. The mixture method involves assigning gRNAs to cells using a principled mixture model. Next, the thresholding method assigns a gRNA to a cell if the UMI count of the gRNA in the cell is greater than or equal to some integer threshold. Finally, the maximum method assigns the gRNA that accounts for the greatest number of UMIs in a given cell to that cell. The maximum method is available only in low MOI. See \href{https://timothy-barry.github.io/sceptre-book/assign-grnas.html}{Chapter 3 of the manual} for more detailed information about \code{assign_grnas()}.
}
\note{
See the manual for information about the method-specific additional arguments.
}
\examples{
library(sceptredata)
data("lowmoi_example_data")
# 1. import data, set default analysis parameters
sceptre_object <- import_data(
  response_matrix = lowmoi_example_data$response_matrix,
  grna_matrix = lowmoi_example_data$grna_matrix,
  extra_covariates = lowmoi_example_data$extra_covariates,
  grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
  moi = "low"
) |> set_analysis_parameters()

# 2. assign gRNAs (three different methods)
sceptre_object <- sceptre_object |> assign_grnas(method = "thresholding")
sceptre_object <- sceptre_object |> assign_grnas(method = "maximum")
sceptre_object <- sceptre_object |> assign_grnas(
  method = "mixture", parallel = TRUE, n_processors = 2
)
}
