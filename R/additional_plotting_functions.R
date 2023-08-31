#' Plot resampling distribution
#'
#' Plots the resampling distribution of null test statistics for a given response-gRNA group pair.
#'
#' The fitted parametric curve is superimposed over the histogram in purple (if applicable). The original test statistic, meanwhile, is plotted as a vertical blue line.
#'
#' @param result_df the output of a call to \code{run_sceptre_lowmoi} (or \code{run_sceptre_highmoi_experimental}). NOTE: the argument \code{output_amount} should be set to \code{3} so that the null statistics are returned.
#' @param row_number an integer identifying the row of the results data frame to visualize
#'
#' @return a plot of the resampling distribution of null test statistics.
#' @export
#'
#' @examples
#' library(Matrix)
#' set.seed(5)
#'
#' \dontrun{
#' # 0. load the data associated with the experiment
#' data(response_matrix_lowmoi) # response-by-cell expression matrix
#' data(grna_matrix_lowmoi) # gRNA-by-cell expression matrix
#' data(covariate_data_frame_lowmoi) # cell-by-covariate data frame
#' data(grna_group_data_frame_lowmoi) # gRNA group information
#' response_grna_group_pairs <- generate_all_pairs(response_matrix_lowmoi,
#' grna_group_data_frame_lowmoi) |> dplyr::sample_n(30)
#' formula_object <- formula(~log(response_n_umis) + log(response_n_nonzero) +
#' bio_rep + p_mito)
#'
#' # 1. apply sceptre
#' discovery_result <- run_sceptre_lowmoi(response_matrix = response_matrix_lowmoi,
#' grna_matrix = grna_matrix_lowmoi,
#' covariate_data_frame = covariate_data_frame_lowmoi,
#' grna_group_data_frame = grna_group_data_frame_lowmoi,
#' formula_object = formula_object,
#' response_grna_group_pairs = response_grna_group_pairs,
#' calibration_check = FALSE,
#' output_amount = 3)
#'
#' # 2. visualize a resampling distribution
#' plot_resampling_distribution(discovery_result, 2)
#' }
plot_resampling_distribution <- function(result_df, row_number) {
  # check for the presence of z_null_1, z_null_2, z_null_3, ...
  df_colnames <- colnames(result_df)
  null_dist_idxs <- grep(pattern = "^z_null", x = df_colnames)
  if (length(null_dist_idxs) == 0) stop("Columns z_null_1, z_null_2, z_null_3, ... must be present in `result_df` to use this function. Rerun your calibration check or discovery analysis, setting `output_amount` to 3.")

  # extract the null and observed z-scores
  z_null <- as.numeric(result_df[row_number, null_dist_idxs]) |> stats::na.omit()
  z_obs <- as.numeric(result_df[row_number, "z_orig"])

  # extract the fitted skew-normal parameters
  fitted_params <- as.numeric(result_df[row_number, c("xi", "omega", "alpha")])
  plot_density <- all(!is.na(fitted_params))

  # create the plot
  p <- data.frame(z_null = z_null) |>
    ggplot2::ggplot(ggplot2::aes(x = z_null)) +
    ggplot2::geom_histogram(ggplot2::aes(y =  ggplot2::after_stat(density)),
                            bins = 50, color = "black", fill = "grey80") +
    ggplot2::geom_vline(xintercept = z_obs, col = "blue", lwd = 0.7) +
    (if (plot_density) ggplot2::stat_function(fun = function(x) {
      compute_sn_density(x, fitted_params[1], fitted_params[2], fitted_params[3])
    }, col = "purple", lwd = 0.7) else NULL) +
    ggplot2::xlab("Null distribution") + ggplot2::ylab("Density") + get_my_theme() +
    ggplot2::scale_y_continuous(expand = c(0, 0))

  return(p)
}
