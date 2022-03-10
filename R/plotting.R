#' Plot fitted skew-t distribution
#'
#' Plots the skew-t distribution fitted to the resampled test statistics, alongside the "ground truth" test statistic derived from the raw data.
#'
#' One also can call this function on a given row of the data frame outputted by `run_sceptre_in_memory` when `full_output` is set to TRUE.
#'
#' @param sceptre_result output of `run_sceptre_gRNA_gene_pair` when `full_output` is set to TRUE
#' @param interval (optional) interval over which the distribution is plotted
#'
#' @return a ggplot object
#' @examples
#' data(gene_matrix); data(gRNA_matrix); data(covariate_matrix)
#' gene_expressions <- gene_matrix[1,]
#' gRNA_expressions <- gRNA_matrix[1,]
#' # run method
#' result <- run_sceptre_gRNA_gene_pair(gene_expressions, gRNA_expressions,
#' covariate_matrix, "left", full_output = TRUE)
#' # plot result
#' plot_skew_t(result)
plot_skew_t <- function(sceptre_result, interval = c(-4, 4)) {
  resampled_zvalues <- sceptre_result %>% dplyr::select(dplyr::starts_with("z_null_")) %>% as.numeric()
  original_zvalue <- sceptre_result$z_value
  dp <- sceptre_result %>% dplyr::select(xi, omega, alpha, nu) %>% as.numeric()
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
