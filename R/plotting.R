#' Plot resampled test statistics over fitted skew-t distribution
#'
#' Plots the fitted skew-t distribution, along with the resampled (null) test statistics and real test statistic.
#'
#' @param resampled_zvalues the set of resampled z-values
#' @param original_zvalue the z value of the original negative binomial fit
#' @param dp the skew-t fit MLE
#' @param interval (optional; default c(-4,4)) interval over which to plot the skew-t distribution
#'
#' @return a ggplot object containing the plot
#' @export
#' @examples
#' data(expressions)
#' data(gRNA_indicators)
#' data(covariate_matrix)
#' result <- run_sceptre_gRNA_gene_pair(expressions = expressions,
#' gRNA_indicators = gRNA_indicators,
#' covariate_matrix = covariate_matrix,
#' reduced_output = FALSE,
#' seed = 4,
#' verbose = FALSE)
#' plot_skew_t(result$resampled_z_values, result$z_value, result$skew_t_mle)
plot_skew_t <- function(resampled_zvalues, original_zvalue, dp, interval = c(-4,4)) {
  z <- seq(interval[1], interval[2], length.out = 1000)
  df_curves <- data.frame(z = z, fitted = sn::dst(x = z, dp = dp), gaussian = stats::dnorm(z)) %>%
    tidyr::gather("curve", "y", fitted, gaussian) %>%
    dplyr::mutate(curve = factor(curve, levels = c("fitted","gaussian"), labels = c("Conditional\nrandomization","Negative\nbinomial")))
  df_ribbon <- df_curves %>%
    dplyr::filter(z <= original_zvalue, curve == "Conditional\nrandomization") %>%
    dplyr::mutate(lower = 0, upper = y) %>%
    dplyr::select(z, lower, upper)
  p <- ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(x = z, y = ..density..),
                   data = data.frame(z = resampled_zvalues),
                   boundary = 0, colour = "black", fill = "lightgray", bins = 15) +
    ggplot2::geom_line(ggplot2::aes(x = z, y = y, group = curve, colour = curve, linetype = curve), data = df_curves) +
    ggplot2::geom_vline(xintercept = original_zvalue, colour = "firebrick3", linetype = "solid") +
    ggplot2::geom_ribbon(ggplot2::aes(x = z, ymin = lower, ymax = upper), fill = "darkorchid2", alpha = 0.5, data = df_ribbon) +
    ggplot2::scale_colour_manual(values = c("darkorchid2", "black"), name = "Null distribution") +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed"), name = "Null distribution") +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::xlab("Negative binomial z-value") + ggplot2::theme_bw() +
    ggplot2::theme(legend.position = c(0.85, 0.8),
          legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
          plot.title = ggplot2::element_text(hjust = 0.5, size = 11),
          panel.grid = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_line(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank())
  return(p)
}
