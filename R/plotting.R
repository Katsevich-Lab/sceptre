#' Plot result
#'
#' For a given gRNA-gene pair analyzed by `sceptre`, plots the resampled test statistics alongside the "ground truth" test statistic derived from the raw data.
#'
#' @param row single row of the data frame outputted by `run_sceptre_high_moi`, when `full_output` is set to TRUE
#'
#' @return a ggplot2 object containing the plot
#' @export
#' @examples
#' # RUN THE METHOD
#' set.seed(4)
#' library(magrittr)
#' library(dplyr)
#' data(gene_matrix)
#' data(gRNA_matrix)
#' data(covariate_matrix)
#' data(gRNA_groups_table)
#' data(gene_gRNA_group_pairs)
#' combined_perturbation_matrix <- threshold_gRNA_matrix(gRNA_matrix) %>%
#' combine_perturbations(gRNA_groups_table)
#' gene_gRNA_group_pairs <- gene_gRNA_group_pairs %>% sample_n(1)
#' result <- run_sceptre_high_moi(gene_matrix = gene_matrix,
#' combined_perturbation_matrix = combined_perturbation_matrix,
#' covariate_matrix = covariate_matrix,
#' gene_gRNA_group_pairs = gene_gRNA_group_pairs,
#' side = "left", parallel = FALSE, full_output = TRUE)
#'
#' # CREATE THE PLOT
#' plot_result(result[1,])
plot_result <- function(row) {
  resampled_zvalues <- row %>% dplyr::select(dplyr::starts_with("z_null_")) %>% as.numeric()
  original_zvalue <- row$z_value
  interval <- range(c(resampled_zvalues, original_zvalue)) + c(-0.5, 0.5)
  dp <- row %>% dplyr::select(xi, omega, alpha, nu) %>% as.numeric()
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


#' Make QQ-plot
#'
#' Makes a QQ-plot for a set of values hypothesized to follow a uniform distribution (e.g., *p*-values).
#'
#' @param p_values a vector values -- hypothesized to follow a uniform distribution -- from which to construct the QQ-plot. This vector typically will be a set of *p*-values.
#' @param ci_level level of the pointwise confidence band (default 0.95)
#' @param point_col color of the plotted points
#' @param alpha transparency of the plotted points
#'
#' @return a ggplot object containing the QQ-plot
#' @export
#' @examples
#' set.seed(4)
#' p_vals <- runif(5000)
#' make_qq_plot(p_vals)
make_qq_plot <- function(p_values, ci_level = 0.95, point_col = "royalblue4", alpha = 0.9) {
  p_thresh <- 1e-8
  to_plot <- data.frame(pvalue = p_values) %>%
    dplyr::mutate(r = rank(pvalue), expected = stats::ppoints(dplyr::n())[r],
                  clower = stats::qbeta(p = (1 - ci_level)/2, shape1 = r, shape2 = dplyr::n() + 1 - r),
                  cupper = stats::qbeta(p = (1 + ci_level)/2, shape1 = r, shape2 = dplyr::n()+ 1 - r)) %>%
    # dplyr::filter(-log10(expected) > 2) %>%
    dplyr::mutate(pvalue = ifelse(pvalue < p_thresh, p_thresh, pvalue))
  p <- ggplot2::ggplot(data = to_plot, mapping = ggplot2::aes(x = expected, y = pvalue, ymin = clower, ymax = cupper)) +
    ggplot2::geom_point(size = 1, alpha = alpha, col = point_col) +
    ggplot2::geom_ribbon(alpha = 0.25) +
    ggplot2::geom_abline(intercept = 0, slope = 1) +
    ggplot2::scale_x_continuous(trans = revlog_trans(base = 10)) + ggplot2::scale_y_continuous(trans = revlog_trans(base = 10)) +
    ggplot2::xlab(expression(paste("Expected null p-value"))) +
    ggplot2::ylab(expression(paste("Observed p-value"))) +
    ggplot2::ggtitle("QQ-plot of p-values") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = c(0.25, 0.85),
                       legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
                       legend.title = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank(),
                       strip.background = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       axis.line = ggplot2::element_line(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::geom_vline(xintercept = 0.1, col = "darkred", linetype = "dashed")
  return(p)
}


StatQQBand <- ggplot2::ggproto("StatQQBand", ggplot2::Stat,
  required_aes = c("y"),

  setup_data = function(data, params){
    if ("group" %in% names(data)) {
      max_unique_pts_per_group <- data |>
        dplyr::group_by(PANEL, group) |>
        dplyr::summarise(pts_per_group = dplyr::n()) |>
        dplyr::summarise(unique_pts_per_group = length(unique(pts_per_group))) |>
        dplyr::summarise(max(unique_pts_per_group)) |>
        dplyr::pull()
      if (max_unique_pts_per_group > 1) {
        stop("Within each panel, you must have the same number of points per QQ plot!")
      } else {
        data |>
          dplyr::select(y, PANEL, group) |>   # remove attributes like color, shape, etc.
          dplyr::filter(group == min(group))  # keep only one of the groups to plot
      }
    } else {
      data
    }
  },

  compute_group = function(data, scales, distribution = "unif", ci_level = 0.95) {
    quantile_fun <- eval(parse(text = sprintf("stats::q%s", distribution)))
    if(is.null(scales$x$trans)){
      scales$x$trans <- scales::identity_trans()
    }
    if(is.null(scales$y$trans)){
      scales$y$trans <- scales::identity_trans()
    }
    data |>
      dplyr::mutate(y = scales$y$trans$inverse(y)) |>
      dplyr::mutate(r = rank(y),
                    x = quantile_fun(stats::ppoints(dplyr::n())[r]),
                    ymin = quantile_fun(stats::qbeta(p = (1 - ci_level)/2, shape1 = r, shape2 = dplyr::n() + 1 - r)),
                    ymax = quantile_fun(stats::qbeta(p = (1 + ci_level)/2, shape1 = r, shape2 = dplyr::n()+ 1 - r))) |>
      dplyr::mutate(x = scales$x$trans$transform(x),
                    y = scales$y$trans$transform(y),
                    ymin = scales$y$trans$transform(ymin),
                    ymax = scales$y$trans$transform(ymax))
  }
)

stat_qq_band <- function(mapping = NULL, data = NULL, geom = "ribbon",
                           position = "identity", na.rm = FALSE, show.legend = FALSE,
                           inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatQQBand, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, alpha = 0.25, ...)
  )
}

StatQQPoints <- ggplot2::ggproto("StatQQPoints", ggplot2::Stat,
  required_aes = c("y"),

  compute_group = function(data, scales, distribution = "unif") {
    if(is.null(scales$x$trans)){
      scales$x$trans <- scales::identity_trans()
    }
    if(is.null(scales$y$trans)){
      scales$y$trans <- scales::identity_trans()
    }
    quantile_fun <- eval(parse(text = sprintf("stats::q%s", distribution)))
    data |>
      dplyr::mutate(y = scales$y$trans$inverse(y)) |>
      dplyr::mutate(x = quantile_fun(stats::ppoints(dplyr::n())[rank(y)])) |>
      dplyr::mutate(x = scales$x$trans$transform(x),
                    y = scales$y$trans$transform(y))
  }
)

stat_qq_points <- function(mapping = NULL, data = NULL, geom = "point",
                           position = "identity", na.rm = FALSE, show.legend = NA,
                           inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatQQPoints, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

data <- tidyr::expand_grid(tibble::tibble(dataset_A = c(rep("A1", 100),
                                                rep("A1", 200),
                                                rep("A2", 300),
                                                rep("A2", 400)),
                                  dataset_B = c(rep("B1", 100),
                                                rep("B2", 200),
                                                rep("B1", 300),
                                                rep("B2", 400))),
                   tibble::tibble(method = c("X", "Y"))) |>
  dplyr::mutate(p_value = rnorm(dplyr::n()))

# data <- tibble::tibble(p_vals = runif(200), method = factor(rep(c(1,2), 100)), type = factor(c(rep(1, 100), rep(2, 100))))
data |>
  ggplot2::ggplot(ggplot2::aes(y = p_value, colour = method)) +
  stat_qq_points(distribution = "norm") +
  stat_qq_band(distribution = "norm") +
  ggplot2::geom_abline() +
  # ggplot2::scale_x_continuous(trans = revlog_trans(base = 10)) +
  # ggplot2::scale_y_continuous(trans = revlog_trans(base = 10)) +
  ggplot2::facet_grid(dataset_A ~ dataset_B) +
  ggplot2::theme_bw()

ci_level = 0.95
plot_2 <- tibble::tibble(y = p_vals) |>
  dplyr::mutate(r = rank(y),
                x = stats::ppoints(dplyr::n())[r],
                ymin = stats::qbeta(p = (1 - ci_level)/2, shape1 = r, shape2 = dplyr::n() + 1 - r),
                ymax = stats::qbeta(p = (1 + ci_level)/2, shape1 = r, shape2 = dplyr::n()+ 1 - r)) |>
  ggplot2::ggplot(ggplot2::aes(x = x, y = y, ymin = ymin, ymax = ymax)) +
  ggplot2::geom_point() +
  ggplot2::geom_ribbon(alpha = 0.25) +
  ggplot2::geom_abline() +
  ggplot2::scale_x_continuous(trans = revlog_trans(base = 10)) +
  ggplot2::scale_y_continuous(trans = revlog_trans(base = 10)) +
  ggplot2::theme_bw()
plot(plot_2)

revlog_trans <- function(base = exp(1)) {
  trans <- function(x) {
    -log(x, base)
  }
  inv <- function(x){
    base^(-x)
  }
  scales::trans_new(paste("revlog-", base, sep = ""),
                    trans,
                    inv,
                    scales::log_breaks(base = base),
                    domain = c(1e-100, Inf)
  )
}
