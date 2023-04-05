#' Plot calibration result
#'
#' The \code{plot_calibration_result()} function helps to visualize the results of a calibration check.
#'
#' The plot contains four panels.
#' \itemize{
#' \item{upper left:} a Q-Q plot of the p-values on an untransformed scale. The p-values should lie along the diagonal line.
#' \item{upper right:} a Q-Q plot of the p-values on a negative log-10 transformed scale. The p-values should lie along the diagonal line, and the majority of the p-values should fall within the gray confidence band.
#' \item{bottom left:} a histogram of the estimated log-2 fold changes in expression. The histogram should be roughly symmetric and centered at zero.
#' \item{bottom right:} a text box displaying (i) the number of false discoveries made on the negative control pairs and (ii) the mean estimated log-fold change. The number of false discoveries ideally should be zero, although a small, positive integer (such as 1, 2, or 3) also is OK. The mean estimated log-fold change should be a numeric value close to zero.
#' }
#'
#' @inheritParams compare_calibration_and_discovery_results
#' @param return_indiv_plots (optional; default \code{FALSE}) return a single combined plot (\code{TRUE}) or the individual plots in a list (\code{FALSE})?
#'
#' @return a single \code{cowplot} object containing the combined panels (if \code{return_indiv_plots} is set to \code{TRUE}) or a list of the individual panels (if \code{return_indiv_plots} is set to \code{FALSE}).
#' @export
#'
#' @note
#' \itemize{
#' \item{The untransformed and transformed Q-Q plots are both informative: the former enables us to see the bulk of the p-value distribution, while the latter enables visualize the tail.}
#' \item{The number of false discoveries depends both on \code{alpha} and \code{multiple_testing_correction}. The default value of these arguments is 0.1 and "BH". This corresponds to a BH correction at nominal false discovery rate (FDR) 0.1.}
#' \item{Technical point: when applying BH at level 0.1 to a collection of strictly null p-values, BH controls family-wise error rate (FWER) at level 0.1 as well as FDR at level 0.1. FWER is the probability of making one or more false discoveries. Thus, with probability 0.9, the number of rejections that BH makes on (well-calibrated) null p-values at level 0.1 is 0. This implies that \code{sceptre} (or any method for that matter) should make about zero false discoveries on negative control p-values data after applying a BH correction at level 0.1.}
#' }
#' @examples
#' # See the example in the run_sceptre_lowmoi help file.
#' ?run_sceptre_lowmoi
plot_calibration_result <- function(calibration_result, alpha = 0.1, multiple_testing_correction = "BH", return_indiv_plots = FALSE) {
  # set the theme
  my_theme <- get_my_theme()

  # compute n rejections
  n_rejections <- obtain_discovery_set(calibration_result, alpha, multiple_testing_correction) |> nrow()
  str1 <- paste0("Number of\nfalse discoveries\n(at alpha ", signif(alpha, 1), "): ", n_rejections)
  str2 <- paste0("\n\nMean log-fold\nchange: ", signif(mean(calibration_result$log_2_fold_change), 2))
  str <- paste0(str1, str2)

  p_a <- ggplot2::ggplot(data = calibration_result,
                         mapping = ggplot2::aes(y = p_value)) +
    stat_qq_points(ymin = 1e-8, size = 0.55,
                                col = "firebrick2",
                                alpha = 0.9) +
    stat_qq_band() +
    ggplot2::scale_x_reverse() +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(x = "Expected null p-value", y = "Observed p-value") +
    ggplot2::geom_abline(col = "black") +
    ggplot2::ggtitle("Q-Q plot (untransformed)") +
    my_theme

  p_b <- ggplot2::ggplot(data = calibration_result, mapping = ggplot2::aes(y = p_value)) +
    stat_qq_points(ymin = 1e-8, size = 0.55,
                   col = "firebrick2", alpha = 0.9) +
    stat_qq_band() +
    ggplot2::scale_x_continuous(trans = revlog_trans(10)) +
    ggplot2::scale_y_continuous(trans = revlog_trans(10)) +
    ggplot2::labs(x = "Expected null p-value", y = "Observed p-value") +
    ggplot2::geom_abline(col = "black") +
    ggplot2::ggtitle("Q-Q plot (transformed)") +
    my_theme

  p_c <- ggplot2::ggplot(data = calibration_result |> dplyr::filter(abs(log_2_fold_change) < 0.6),
                mapping = ggplot2::aes(x = log_2_fold_change, ggplot2::after_stat(density))) +
    ggplot2::geom_histogram(binwidth = 0.05, fill = "grey90", col = "black", boundary = 0) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.0, .01))) +
    ggplot2::ggtitle("Log fold changes") +
    ggplot2::xlab("Estimated log fold change") + ggplot2::ylab("Density") +
    ggplot2::geom_vline(xintercept = 0, col = "firebrick2", linewidth = 1) +
    my_theme

  p_d <- ggplot2::ggplot() +
    ggplot2::annotate(geom = "text", label = str, x = 1.3, y = 1.2) +
    ggplot2::theme_void() +
    ggplot2::xlim(c(0, 2)) +
    ggplot2::ylim(c(0, 2))

  if (return_indiv_plots) {
    out <- list(p_a, p_b, p_c, p_d)
  } else {
    out <- cowplot::plot_grid(p_a, p_b, p_c, p_d, nrow = 2, rel_heights = c(0.55, 0.45))
  }
  return(out)
}


#' Compare calibration and discovery results
#'
#' The \code{compare_calibration_and_discovery_results()} function compares two sets of p-values: the negative control p-values (i.e., those obtained from the calibration check) and the discovery p-values (i.e., those obtained from the discovery analysis). These two sets of p-values are plotted on the same Q-Q plot, with the negative control p-values colored in red and the discovery p-values colored in blue. The negative control p-values should lie along the diagonal line and fall mostly within the gray confidence band, indicating control of the false discovery rate. By contrast, the discovery p-values should lie above the diagonal line, indicating the presence of signal in the discovery set.
#'
#' The horizontal dashed line indicates the multiple testing threshold. Blue points above the horizontal dashed line are called as discoveries.
#'
#' @param calibration_result (required) the output of a call to \code{run_sceptre_lowmoi} with \code{calibration_check} set to \code{TRUE}
#' @param discovery_result (required) the output of a call to \code{run_sceptre_lowmoi} with \code{calibration_check} set to \code{FALSE}
#' @param alpha (optional; default 0.1) the nominal type-I error control level
#' @param multiple_testing_correction (optional; default "BH") the multiple testing correction method (one of "BH", "bonferroni", "BY", or "holm")
#'
#' @return A ggplot2 object containing the Q-Q plot.
#' @export
#'
#' @examples
#' # See the example in the run_sceptre_lowmoi help file.
#' ?run_sceptre_lowmoi
compare_calibration_and_discovery_results <- function(calibration_result, discovery_result, alpha = 0.1, multiple_testing_correction = "BH") {
  discovery_result <- discovery_result |> stats::na.omit()
  lab <- c(rep(factor("Negative control"), nrow(calibration_result)),
           rep(factor("Discovery"), nrow(discovery_result))) |>
    factor(levels = c("Discovery", "Negative control"))

  discovery_set <- obtain_discovery_set(discovery_result, alpha, multiple_testing_correction)
  thresh <- if (nrow(discovery_set) >= 1) max(discovery_set$p_value) else NULL

  df <- data.frame(p_value = c(calibration_result$p_value, discovery_result$p_value),
                   lab = lab)
  ggplot2::ggplot(data = df, mapping = ggplot2::aes(y = p_value, col = lab)) +
    stat_qq_points(ymin = 1e-8, size = 0.8,
                                alpha = 0.8) +
    stat_qq_band() +
    ggplot2::scale_x_continuous(trans = revlog_trans(10)) +
    ggplot2::scale_y_continuous(trans = revlog_trans(10)) +
    ggplot2::labs(x = "Expected null p-value", y = "Observed p-value") +
    ggplot2::geom_abline(col = "black") +
    ggplot2::ggtitle("Q-Q plot (transformed)") +
    get_my_theme(12) +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = "bottom") +
    ggplot2::geom_hline(yintercept = thresh, linetype = "dashed") +
    ggplot2::scale_color_manual(values = c("dodgerblue3", "firebrick2"))
}


get_my_theme <- function(element_text_size = 11) {
  ggplot2::theme_bw() + ggplot2::theme(axis.line = ggplot2::element_line(color = "black"),
                                       panel.grid.major = ggplot2::element_blank(),
                                       panel.grid.minor = ggplot2::element_blank(),
                                       panel.border = ggplot2::element_blank(),
                                       panel.background = ggplot2::element_blank(),
                                       plot.title = ggplot2::element_text(hjust = 0.5, size = element_text_size),
                                       legend.text = ggplot2::element_text(size = 11))
}


#' Make volcano plot
#'
#' The \code{make_volcano_plot()} function creates a volcano plot of the discovery results. Each point in the plot corresponds to a response-gRNA group pair; the estimated log-2 fold change of the pair is plotted on the x-axis, and the (negative log-10 transformed) p-value is plotted on the y-axis.
#'
#'  The horizontal dashed line indicates the multiple testing threshold. Points above the dashed line (which are colored in purple) are called as discoveries; points below the dashed line (which are colored in blue), meanwhile, are called as insignificant.
#'
#' @inheritParams compare_calibration_and_discovery_results
#'
#' @return a ggplot2 object
#' @export
#' @examples
#' # See the example in the run_sceptre_lowmoi help file.
#' ?run_sceptre_lowmoi
make_volcano_plot <- function(discovery_result, alpha = 0.1, x_limits = c(-1.5, 1.5), transparency = 0.5, point_size = 0.9, multiple_testing_correction = "BH") {
  discovery_result <- discovery_result |> stats::na.omit()
  disc_set <- obtain_discovery_set(discovery_result, alpha, multiple_testing_correction)
  p_thresh <- max(disc_set$p_value)
  p_lower_lim <- 1e-12
  out <- ggplot2::ggplot(data = discovery_result |> dplyr::mutate(reject = p_value < p_thresh,
                                                                  p_value = ifelse(p_value < p_lower_lim, p_lower_lim, p_value),
                                                                  log_2_fold_change = ifelse(log_2_fold_change > x_limits[2], x_limits[2], log_2_fold_change),
                                                                  log_2_fold_change = ifelse(log_2_fold_change < x_limits[1], x_limits[1], log_2_fold_change)),
                  mapping = ggplot2::aes(x = log_2_fold_change, y = p_value, col = reject)) +
    ggplot2::geom_point(alpha = transparency, size = point_size) +
    ggplot2::scale_y_continuous(trans = revlog_trans(10), limits = c(1, 1e-12), expand = c(0.02, 0)) +
    get_my_theme() + ggplot2::xlab("Log fold change") + ggplot2::ylab("P-value") +
    ggplot2::geom_hline(yintercept = p_thresh, linetype = "dashed") +
    ggplot2::theme(legend.position = "none") + ggplot2::scale_color_manual(values = c("dodgerblue", "blueviolet")) +
    ggplot2::ggtitle("P-value vs. log fold change")
  return(out)
}
