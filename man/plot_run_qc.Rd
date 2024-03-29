% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plot_run_qc}
\alias{plot_run_qc}
\title{Plot run QC}
\usage{
plot_run_qc(
  sceptre_object,
  downsample_pairs = 10000L,
  point_size = 0.55,
  transparency = 0.8,
  return_indiv_plots = FALSE
)
}
\arguments{
\item{sceptre_object}{a \code{sceptre_object} that has had \code{run_qc()} called on it}

\item{downsample_pairs}{(optional; default \code{10000}) the maximum number of points to plot in the lower panel of the figure (i.e., the "pairwise QC" plot)}

\item{point_size}{(optional; default \code{0.55}) the size of the individual points in the plot}

\item{transparency}{(optional; default \code{0.8}) the transparency of the individual points in the plot}

\item{return_indiv_plots}{(optional; default \code{FALSE}) if \code{FALSE} then a list of \code{ggplot} is returned; if \code{TRUE} then a single \code{cowplot} object is returned.}
}
\value{
a single \code{cowplot} object containing the combined panels (if \code{return_indiv_plots} is set to \code{TRUE}) or a list of the individual panels (if \code{return_indiv_plots} is set to \code{FALSE})
}
\description{
\code{plot_run_qc()} creates a visualization of the outcome of the QC step. See \href{https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_qc}{Section 4 of the introductory chapter in the manual} for guidance on interpreting this visual.
}
\examples{
# A full example can be found at ?sceptre;
# `plot_run_qc()` is dispatched when
# `plot()` is called on the `sceptre_object`
# in step 4 (the run qc step).
}
