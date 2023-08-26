get_my_theme <- function(element_text_size = 11) {
  ggplot2::theme_bw() + ggplot2::theme(axis.line = ggplot2::element_line(color = "black"),
                                       panel.grid.major = ggplot2::element_blank(),
                                       panel.grid.minor = ggplot2::element_blank(),
                                       panel.border = ggplot2::element_blank(),
                                       panel.background = ggplot2::element_blank(),
                                       plot.title = ggplot2::element_text(hjust = 0.5, size = element_text_size))
}

#' Plot calibration result
#'
#' The \code{plot_calibration_result()} function helps to visualize the results of a calibration check.
#'
#' The plot contains four panels.
#' \itemize{
#' \item{upper left:} a QQ plot of the p-values on an untransformed scale. The p-values should lie along the diagonal line.
#' \item{upper right:} a QQ plot of the p-values on a negative log-10 transformed scale. The p-values should lie along the diagonal line, and the majority of the p-values should fall within the gray confidence band.
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
#' \item{The untransformed and transformed QQ plots are both informative: the former enables us to see the bulk of the p-value distribution, while the latter enables visualize the tail.}
#' \item{The number of false discoveries depends both on \code{alpha} and \code{multiple_testing_correction}. The default value of these arguments is 0.1 and "BH". This corresponds to a BH correction at nominal false discovery rate (FDR) 0.1.}
#' \item{Technical point: when applying BH at level 0.1 to a collection of strictly null p-values, BH controls family-wise error rate (FWER) at level 0.1 as well as FDR at level 0.1. FWER is the probability of making one or more false discoveries. Thus, with probability 0.9, the number of rejections that BH makes on (well-calibrated) null p-values at level 0.1 is 0. This implies that \code{sceptre} (or any method for that matter) should make about zero false discoveries on negative control p-values data after applying a BH correction at level 0.1.}
#' }
#' @examples
#' # See the example in the run_sceptre_lowmoi help file.
#' ?run_sceptre_lowmoi
plot_calibration_result <- function(sceptre_object, return_indiv_plots = FALSE, point_size = 0.55, transparency = 0.8) {
  calibration_result <- sceptre_object@calibration_result
  if (nrow(calibration_result) == 0L) stop("Calibration check not yet called.")
  my_theme <- get_my_theme()

  # compute n rejections
  n_rejections <- sum(calibration_result$reject)
  str1 <- paste0("Number of\nfalse discoveries\n(at alpha ", signif(sceptre_object@multiple_testing_alpha, 1), "): ", n_rejections)
  str2 <- paste0("\n\nMean log-fold\nchange: ", signif(mean(calibration_result$log_2_fold_change), 2))
  str <- paste0(str1, str2)

  p_a <- ggplot2::ggplot(data = calibration_result,
                         mapping = ggplot2::aes(y = p_value)) +
    stat_qq_points(ymin = 1e-8, size = point_size,
                                col = "firebrick2",
                                alpha = transparency) +
    stat_qq_band() +
    ggplot2::scale_x_reverse() +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(x = "Expected null p-value", y = "Observed p-value") +
    ggplot2::geom_abline(col = "black") +
    ggplot2::ggtitle("QQ plot (bulk)") +
    my_theme

  p_b <- ggplot2::ggplot(data = calibration_result, mapping = ggplot2::aes(y = p_value)) +
    stat_qq_points(ymin = 1e-8, size = point_size,
                   col = "firebrick2", alpha = transparency) +
    stat_qq_band() +
    ggplot2::scale_x_continuous(trans = revlog_trans(10)) +
    ggplot2::scale_y_continuous(trans = revlog_trans(10)) +
    ggplot2::labs(x = "Expected null p-value", y = "Observed p-value") +
    ggplot2::geom_abline(col = "black") +
    ggplot2::ggtitle("QQ plot (tail)") +
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

compare_calibration_and_discovery_results <- function(calibration_result, discovery_result, p_thresh,
                                                      transform_scale = TRUE, include_legend = FALSE,
                                                      include_y_axis_text = TRUE, point_size = 0.55,
                                                      transparency = 0.8) {
  lab <- c(rep(factor("Negative control"), nrow(calibration_result)),
           rep(factor("Discovery"), nrow(discovery_result))) |>
    factor(levels = c("Discovery", "Negative control"))
  df <- data.frame(p_value = c(calibration_result$p_value, discovery_result$p_value), lab = lab)

  p_out <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(y = p_value, col = lab)) +
    stat_qq_points(ymin = 1e-8, size = point_size, alpha = transparency) +
    stat_qq_band() +
    ggplot2::labs(x = "Expected null p-value", y = "Observed p-value") +
    ggplot2::geom_abline(col = "black") +
    get_my_theme() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = if (include_legend) "bottom" else "none",
                   axis.title.y = if (include_y_axis_text) ggplot2::element_text() else ggplot2::element_blank()) +
    ggplot2::scale_color_manual(values = c("dodgerblue3", "firebrick2"))

  if (!transform_scale) {
    p_out <- p_out +
      ggplot2::scale_x_reverse() +
      ggplot2::scale_y_reverse() +
      ggplot2::ggtitle("QQ plot (bulk)")
  } else {
    p_out <- p_out +
      ggplot2::scale_x_continuous(trans = revlog_trans(10)) +
      ggplot2::scale_y_continuous(trans = revlog_trans(10)) +
      ggplot2::geom_hline(yintercept = p_thresh, linetype = "dashed") +
      ggplot2::ggtitle("QQ plot (tail)")
  }

  if (include_legend) {
    p_out <- p_out +
    ggplot2::theme(legend.position = c(0.7, 0.1),
                   legend.margin = ggplot2::margin(t = -0.5, unit = "cm"),
                   legend.title = ggplot2::element_blank()) +
      ggplot2::guides(color = ggplot2::guide_legend(
        keywidth = 0.0,
        keyheight = 0.15,
        default.unit = "inch",
        override.aes = list(size = 1.25)))
  }

  return(p_out)
}

make_volcano_plot <- function(discovery_result, p_thresh, x_limits = c(-1.5, 1.5), transparency = 0.5, point_size = 0.55) {
  p_lower_lim <- 1e-12
  out <- ggplot2::ggplot(data = discovery_result |> dplyr::mutate(reject = p_value < p_thresh,
                                                                  p_value = ifelse(p_value < p_lower_lim, p_lower_lim, p_value),
                                                                  log_2_fold_change = ifelse(log_2_fold_change > x_limits[2], x_limits[2], log_2_fold_change),
                                                                  log_2_fold_change = ifelse(log_2_fold_change < x_limits[1], x_limits[1], log_2_fold_change)),
                  mapping = ggplot2::aes(x = log_2_fold_change, y = p_value, col = reject)) +
    ggplot2::geom_point(alpha = transparency, size = point_size) +
    ggplot2::scale_y_continuous(trans = revlog_trans(10), limits = c(1, 1e-12), expand = c(0.02, 0)) +
    get_my_theme() + ggplot2::xlab("Log fold change") + ggplot2::ylab("P-value") +
    (if (!is.na(p_thresh)) ggplot2::geom_hline(yintercept = p_thresh, linetype = "dashed") else NULL) +
    ggplot2::theme(legend.position = "none") + ggplot2::scale_color_manual(values = c("dodgerblue", "blueviolet")) +
    ggplot2::ggtitle("Discovery volcano plot")
  return(out)
}


plot_discovery_result <- function(sceptre_object, return_indiv_plots = FALSE, x_limits = c(-1.5, 1.5), transparency = 0.8, point_size = 0.55) {
  # first, compute the rejection set
  if (nrow(sceptre_object@discovery_result) == 0L) stop("Discovery analysis not run.")
  discovery_result <- sceptre_object@discovery_result |> na.omit()
  calibration_result <- sceptre_object@calibration_result
  discovery_set <- discovery_result |> dplyr::filter(reject)
  p_thresh <- if (nrow(discovery_set) >= 1L) max(discovery_set$p_value) else NA
  # create the bulk QQ plot
  p1 <- compare_calibration_and_discovery_results(calibration_result = calibration_result,
                                                  discovery_result = discovery_result,
                                                  p_thresh = p_thresh,
                                                  point_size = point_size,
                                                  transparency = transparency,
                                                  transform_scale = FALSE,
                                                  include_legend = TRUE,
                                                  include_y_axis_text = TRUE)
  # create tail qq plot
  p2 <- compare_calibration_and_discovery_results(calibration_result = calibration_result,
                                                  discovery_result = discovery_result,
                                                  p_thresh = p_thresh,
                                                  point_size = point_size,
                                                  transparency = transparency,
                                                  transform_scale = TRUE,
                                                  include_legend = FALSE,
                                                  include_y_axis_text = FALSE)
  # make the volcano plot
  p3 <- make_volcano_plot(discovery_result = discovery_result,
                          p_thresh = p_thresh,
                          transparency = transparency,
                          point_size = point_size,
                          x_limits = x_limits)
  # create the final panel
  n_rejections <- nrow(discovery_set)
  n_pairs <- nrow(discovery_result)
  str <- paste0("Number of discovery pairs\nrejected (at alpha ", signif(sceptre_object@multiple_testing_alpha, 1), "):\n", n_rejections, " of ", n_pairs)
  p4 <- ggplot2::ggplot() +
    ggplot2::annotate(geom = "text", label = str, x = 1.1, y = 1.2) +
    ggplot2::theme_void() +
    ggplot2::xlim(c(0, 2)) +
    ggplot2::ylim(c(0, 2))

  # combine the plots
  if (return_indiv_plots) {
    p_out <- list(p_a, p_b, p_c, p_d)
  } else {
    p_out <- cowplot::plot_grid(p1, p2, p3, p4, labels = c("a", "b", "c"), rel_heights = c(0.55, 0.45), nrow = 2)
  }
  return(p_out)
}
