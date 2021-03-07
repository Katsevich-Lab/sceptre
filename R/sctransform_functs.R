#### NOTE: These functions are taken from the sctransform package. Credit goes to Christoph Hafemeister and Rahul Satija. We had to rework these functions slightly to improve their scalability.

#' Regularize thetas
#'
#' @param genes_log_gmean the log-transformed vector of the geometric mean of the gene expressions
#' @param theta the vector of unregularized gene sizes
#' @param theta_regularization transformation to apply to the thetas before fitting nonparametric regression; options "log_theta" and "od_factor;" defaults to "log_theta"
#' @param bw_adjust amount of regularization (greater value -> more regularization). 0 corresponds to no regularization at all. Default value 3.
#' @param plot_me create informative plots? (default FALSE)
#'
#' @return the regularized vector of thetas
#' @noRd
#'
#' @examples
#' genes_log_gmean <- c(runif(1000, 2, 5))
#' names(genes_log_gmean) <- paste0("gene_", 1:1000)
#' theta <- genes_log_gmean * 0.5 + rnorm(100)
#' theta_regularization <- "log_theta"
#' bw_adjust <- 3
#' theta_reg <- regularize_thetas(genes_log_gmean, theta, theta_regularization, 0, TRUE)
regularize_thetas <- function(genes_log_gmean, theta, theta_regularization = "log_theta", bw_adjust = 3, plot_me = FALSE) {
  if (bw_adjust == 0) {
    theta_out <- theta
  } else {
    theta <- pmax(theta, 1e-7)
    dispersion_par <- switch(theta_regularization, log_theta = log10(theta), od_factor = log10(1 + 10^genes_log_gmean/theta), stop("theta_regularization ", theta_regularization, " unknown - only log_theta and od_factor supported at the moment"))
    outliers <- is_outlier(dispersion_par, genes_log_gmean)

    dispersion_par <- dispersion_par[!outliers]
    genes_log_gmean_step1 <- genes_log_gmean[!outliers]

    bw <- stats::bw.SJ(genes_log_gmean_step1) * bw_adjust
    x_points <- pmax(genes_log_gmean, min(genes_log_gmean_step1))
    x_points <- pmin(x_points, max(genes_log_gmean_step1))
    o <- order(x_points)
    fitted_dispersion_par <- numeric(length = length(x_points))
    fitted_dispersion_par[o] <- stats::ksmooth(x = genes_log_gmean_step1, y = dispersion_par, x.points = x_points, bandwidth = bw, kernel = "normal")$y
    theta_out <- switch(theta_regularization, log_theta = 10^fitted_dispersion_par, od_factor = 10^genes_log_gmean/(10^fitted_dispersion_par - 1))
    problem_idxs <- is.na(theta_out)
    theta_out[problem_idxs] <- theta[problem_idxs]
    attr(theta_out, "outliers") <- outliers
    out <- theta_out

    print(plot_me)
    if (plot_me) {
      df1 <- data.frame(log_theta = c(dispersion_par, fitted_dispersion_par), mean_gene_exp = c(genes_log_gmean_step1, x_points), regularized = c(rep(FALSE, length(genes_log_gmean_step1)), rep(TRUE, length(x_points))))
      p1 <- ggplot2::ggplot(data = df1, mapping = aes(x = mean_gene_exp, y = log_theta, col = regularized)) + geom_point() + theme_bw()
      df2 <- data.frame(theta = theta, theta_out = theta_out)
      p2 <- ggplot2::ggplot(data = df2, mapping = aes(x = theta, y = theta_out)) + geom_point() + theme_bw()
      print(p1); print(p2)
      out <- list(theta_out, p1, p2)
    }
  }
  return(theta_out)
}


#' Internal sctransform function
#' @noRd
is_outlier <- function(y, x, th = 10) {
  bin.width <- (max(x) - min(x)) * stats::bw.SJ(x)/2
  eps <- .Machine$double.eps * 10
  breaks1 <- seq(from = min(x) - eps, to = max(x) + bin.width,
                 by = bin.width)
  breaks2 <- seq(from = min(x) - eps - bin.width/2, to = max(x) +
                   bin.width, by = bin.width)
  score1 <- robust_scale_binned(y, x, breaks1)
  score2 <- robust_scale_binned(y, x, breaks2)
  return(pmin(abs(score1), abs(score2)) > th)
}


#' Internal sctransform function
#' @noRd
robust_scale_binned <- function (y, x, breaks) {
  bins <- cut(x = x, breaks = breaks, ordered_result = TRUE)
  tmp <- aggregate(x = y, by = list(bin = bins), FUN = robust_scale)
  score <- rep(0, length(x))
  o <- order(bins)
  if (inherits(x = tmp$x, what = "list")) {
    score[o] <- unlist(tmp$x)
  }
  else {
    score[o] <- as.numeric(t(tmp$x))
  }
  return(score)
}


#' Internal sctransform function
#' @noRd
robust_scale <- function(x) (x - stats::median(x))/(stats::mad(x) + .Machine$double.eps)


#' Log geometric mean
#'
#' @param x a numeric vector
#' @param eps small number to add to each value (default 1)
#'
#' @return the log-transformed geometric means
#' @noRd
log_geom_mean <- function(x, eps = 1) {log10(exp(mean(log(x + eps))) - eps)}
