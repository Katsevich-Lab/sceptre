# workhorse function for \code{stat_qq_band}
StatQQBand <- ggplot2::ggproto("StatQQBand", ggplot2::Stat,
  required_aes = c("y"), # the sample will come through the y aesthetic
  # code below implements logic that you only want to have one band per panel,
  # rather than having one band per group
  setup_data = function(data, params) {
    if ("group" %in% names(data)) {
      # calculate the maximum across panels of the number of unique numbers of
      # points there are in each group
      max_unique_pts_per_group <- data |>
        dplyr::group_by(PANEL, group) |>
        dplyr::summarise(pts_per_group = dplyr::n()) |>
        dplyr::summarise(unique_pts_per_group = length(unique(pts_per_group))) |>
        dplyr::summarise(max(unique_pts_per_group)) |>
        dplyr::pull()
      if (max_unique_pts_per_group > 1) {
        # throw an error if there is a panel with differing numbers of points per group
        stop("Within each panel, you must have the same number of points per QQ plot!")
      } else {
        data |>
          dplyr::select(y, PANEL, group) |> # remove attributes like color, shape, etc.
          dplyr::filter(group == min(group)) # keep only one of the groups to plot
      }
    } else {
      data
    }
  },
  compute_group = function(data, scales, distribution = "unif", max_pts_to_plot = 500, ci_level = 0.95) {
    # get the quantile function from the stats package
    quantile_fun <- eval(parse(text = sprintf("stats::q%s", distribution)))
    # set x and y transformations to identity if they are not given
    if (is.null(scales$x$trans)) {
      scales$x$trans <- scales::identity_trans()
    }
    if (is.null(scales$y$trans)) {
      scales$y$trans <- scales::identity_trans()
    }
    # compute the upper and lower confidence bands
    data |>
      # the given y is already transformed, so transform it back first
      dplyr::mutate(y = scales$y$trans$inverse(y)) |>
      # apply the formulas for the confidence bands
      dplyr::mutate(
        r = rank(y, ties.method = "first"),
        x = quantile_fun(stats::ppoints(dplyr::n())[r]),
        ymin = quantile_fun(stats::qbeta(p = (1 - ci_level) / 2, shape1 = r, shape2 = dplyr::n() + 1 - r)),
        ymax = quantile_fun(stats::qbeta(p = (1 + ci_level) / 2, shape1 = r, shape2 = dplyr::n() + 1 - r))
      ) |>
      # transform
      dplyr::mutate(
        x = scales$x$trans$transform(x),
        y = scales$y$trans$transform(y),
        ymin = scales$y$trans$transform(ymin),
        ymax = scales$y$trans$transform(ymax)
      ) |> # downsample
      dplyr::mutate(
        x_int = ggplot2::cut_interval(x, n = max_pts_to_plot),
        y_int = ggplot2::cut_interval(y, n = max_pts_to_plot)
      ) |>
      dplyr::group_by(x_int, y_int) |>
      dplyr::slice_sample(n = 1) |>
      dplyr::ungroup()
  }
)

#' Add confidence band for a QQ plot
#'
#' @inheritParams ggplot2::stat_identity
#' @param distribution The distribution with respect to which to make QQ plot
#' (e.g. \code{"unif"}, \code{"norm"}, etc.); default \code{"unif"}
#' @param max_pts_to_plot Points are downsampled for plotting purposes, so that
#' the maximum number of points actually plotted is \code{max_pts_to_plot}.
#' @param ci_level The pointwise confidence level for the QQ plot; default 0.95
#' @noRd
stat_qq_band <- function(mapping = NULL, data = NULL, geom = "ribbon",
                         position = "identity", show.legend = FALSE,
                         inherit.aes = TRUE, distribution = "unif",
                         max_pts_to_plot = 500,
                         ci_level = 0.95, ...) {
  ggplot2::layer(
    stat = StatQQBand, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(
      alpha = 0.25, distribution = distribution,
      max_pts_to_plot = max_pts_to_plot, ci_level = ci_level, ...
    )
  )
}

# workhorse function for \code{stat_qq_points}
StatQQPoints <- ggplot2::ggproto("StatQQPoints", ggplot2::Stat,
  required_aes = c("y"),
  compute_group = function(data, scales, distribution = "unif", max_pts_to_plot = 500, ymin = -Inf, ymax = Inf) {
    # set x and y transformations to identity if they are not given
    if (is.null(scales$x$trans)) {
      scales$x$trans <- scales::identity_trans()
    }
    if (is.null(scales$y$trans)) {
      scales$y$trans <- scales::identity_trans()
    }
    # get the quantile function from the stats package
    quantile_fun <- eval(parse(text = sprintf("stats::q%s", distribution)))
    data |>
      # the given y is already transformed, so transform it back first
      dplyr::mutate(
        y = scales$y$trans$inverse(y),
        # apply the formula
        x = quantile_fun(stats::ppoints(dplyr::n())[rank(y, ties.method = "first")]),
        # threshold
        y = pmin(pmax(y, ymin), ymax),
        # transform
        x = scales$x$trans$transform(x),
        y = scales$y$trans$transform(y)
      ) |>
      # downsample
      dplyr::mutate(
        x_int = ggplot2::cut_interval(x, n = max_pts_to_plot),
        y_int = ggplot2::cut_interval(y, n = max_pts_to_plot)
      ) |>
      dplyr::group_by(x_int, y_int) |>
      dplyr::slice_sample(n = 1) |>
      dplyr::ungroup()
  }
)

#' Add points for a QQ plot
#'
#' @inheritParams ggplot2::stat_identity
#' @param distribution The distribution with respect to which to make QQ plot
#' (e.g. \code{"unif"}, \code{"norm"}, etc.); default \code{"unif"}
#' @param max_pts_to_plot Points are downsampled for plotting purposes, so that
#' the maximum number of points actually plotted is \code{max_pts_to_plot}.
#' @param ymin Samples with value less than \code{ymin} will be truncated at
#' \code{ymin} for plotting.
#' @param ymax Samples with value less than \code{ymax} will be truncated at
#' \code{ymax} for plotting.
#' @noRd
stat_qq_points <- function(mapping = NULL, data = NULL, geom = "point",
                           position = "identity", show.legend = NA,
                           distribution = "unif", max_pts_to_plot = 500,
                           ymin = -Inf, ymax = Inf,
                           inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatQQPoints, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(
      distribution = distribution, max_pts_to_plot = max_pts_to_plot,
      ymin = ymin, ymax = ymax, ...
    )
  )
}


#' Reverse logarithm transformation
#'
#' @param base The base with respect to which to take the logarithm (defaults to e).
#'
#' @return A scale transformation object.
#' @noRd
revlog_trans <- function(base = exp(1)) {
  trans <- function(x) {
    -log(x, base)
  }
  inv <- function(x) {
    base^(-x)
  }
  scales::trans_new(
    name = paste("revlog-", base, sep = ""),
    transform = trans,
    inverse = inv,
    breaks = scales::log_breaks(base = base),
    domain = c(1e-100, Inf)
  )
}
