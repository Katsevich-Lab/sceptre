# workhorse function for \code{stat_qq_band}
StatQQBand <- ggplot2::ggproto(
    "StatQQBand",
    ggplot2::Stat,
    required_aes = c("y"), # the sample will come through the y aesthetic
    # Choose one group's rows for the shared QQ confidence band in each panel
    # data: Plot rows with transformed y values, PANEL and group numbers.
    # params: ggplot statistic settings.
    # Returns: Rows from the lowest-numbered group in each panel, keeping their order.
    #   Groups within a panel must have the same number of points to share a band;
    #   otherwise error. Select separately per panel, since group numbers can differ.
    setup_data = function(data, params) {
        if ("group" %in% names(data)) {
            # Require equal group sizes within each panel before selecting a
            # band.
            max_unique_pts_per_group <- data |>
                dplyr::group_by(PANEL, group) |>
                dplyr::summarise(pts_per_group = dplyr::n()) |>
                dplyr::summarise(
                    unique_pts_per_group = length(unique(pts_per_group))
                ) |>
                dplyr::summarise(max(unique_pts_per_group)) |>
                dplyr::pull()
            if (max_unique_pts_per_group > 1) {
                # throw an error if there is a panel with differing numbers of
                # points per group
                stop(
                    "Within each panel, you must have the same number of points per QQ plot!"
                )
            } else {
                data |>
                    dplyr::select(y, PANEL, group) |> # remove attributes like color, shape, etc.
                    # Each panel has its own group IDs; preserve its first band.
                    dplyr::filter(group == stats::ave(group, PANEL, FUN = min))
            }
        } else {
            data
        }
    },
    # Calculate a pointwise QQ confidence band for one plotted group
    # data,scales: Transformed y values and panel scales. Undo the y
    #   transformation before ranking observations.
    # distribution,max_pts_to_plot,ci_level: Reference distribution, display-grid
    #   resolution and pointwise confidence level. Break rank ties by first occurrence.
    # Returns: x/y/ymin/ymax in plotted coordinates. To reduce overplotting, retain
    #   one sampled row per occupied x/y grid cell, using seed 4 and restoring the
    #   previous random state. max_pts_to_plot controls grid resolution, rather than
    #   imposing a total point limit.
    compute_group = function(
        data,
        scales,
        distribution = "unif",
        max_pts_to_plot = 500,
        ci_level = 0.95
    ) {
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
        withr::with_seed(4, {
            data |>
                # the given y is already transformed, so transform it back first
                dplyr::mutate(y = scales$y$trans$inverse(y)) |>
                # apply the formulas for the confidence bands
                dplyr::mutate(
                    r = rank(y, ties.method = "first"),
                    x = quantile_fun(stats::ppoints(dplyr::n())[r]),
                    ymin = quantile_fun(stats::qbeta(
                        p = (1 - ci_level) / 2,
                        shape1 = r,
                        shape2 = dplyr::n() + 1 - r
                    )),
                    ymax = quantile_fun(stats::qbeta(
                        p = (1 + ci_level) / 2,
                        shape1 = r,
                        shape2 = dplyr::n() + 1 - r
                    ))
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
        })
    }
)

#' Add a pointwise confidence band to a QQ plot
#' @inheritParams ggplot2::stat_identity
#' @param distribution Reference distribution, such as "unif" or "norm"; default "unif".
#' @param max_pts_to_plot Number of display-grid intervals per axis. Keep at most
#'   one point per occupied grid cell; this is not a total point limit.
#' @param ci_level Pointwise confidence level, default 0.95.
#' @return ggplot layer with beta order-statistic bands in each panel's scales.
#'   Break rank ties by first occurrence and restore the previous random state
#'   after sampling points with seed 4.
#' @noRd
stat_qq_band <- function(
    mapping = NULL,
    data = NULL,
    geom = "ribbon",
    position = "identity",
    show.legend = FALSE,
    inherit.aes = TRUE,
    distribution = "unif",
    max_pts_to_plot = 500,
    ci_level = 0.95,
    ...
) {
    ggplot2::layer(
        stat = StatQQBand,
        data = data,
        mapping = mapping,
        geom = geom,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(
            alpha = 0.25,
            distribution = distribution,
            max_pts_to_plot = max_pts_to_plot,
            ci_level = ci_level,
            ...
        )
    )
}

# workhorse function for \code{stat_qq_points}
StatQQPoints <- ggplot2::ggproto(
    "StatQQPoints",
    ggplot2::Stat,
    required_aes = c("y"),
    # Calculate QQ points for one group using its panel's scales
    # data,scales: Transformed y values and panel scales. Undo the y
    #   transformation before ranking, with first-occurrence tie breaking.
    # distribution,max_pts_to_plot: Reference distribution and display-grid resolution.
    # ymin,ymax: Display cutoffs on the original y scale, applied after ranking.
    # Returns: x/y rows in plotted coordinates. Sample one row per occupied grid cell
    #   using seed 4, then restore the caller's previous random state.
    compute_group = function(
        data,
        scales,
        distribution = "unif",
        max_pts_to_plot = 500,
        ymin = -Inf,
        ymax = Inf
    ) {
        # set x and y transformations to identity if they are not given
        if (is.null(scales$x$trans)) {
            scales$x$trans <- scales::identity_trans()
        }
        if (is.null(scales$y$trans)) {
            scales$y$trans <- scales::identity_trans()
        }
        # get the quantile function from the stats package
        quantile_fun <- eval(parse(text = sprintf("stats::q%s", distribution)))
        withr::with_seed(4, {
            data |>
                # the given y is already transformed, so transform it back first
                dplyr::mutate(
                    y = scales$y$trans$inverse(y),
                    # apply the formula
                    x = quantile_fun(stats::ppoints(dplyr::n())[rank(
                        y,
                        ties.method = "first"
                    )]),
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
        })
    }
)

#' Add observed versus expected quantiles to a QQ plot
#' @inheritParams ggplot2::stat_identity
#' @param distribution Reference distribution, such as "unif" or "norm"; default "unif".
#' @param max_pts_to_plot Display-grid intervals per axis; keep at most one point
#'   per occupied grid cell, rather than limiting the total number of points.
#' @param ymin,ymax Original-scale y bounds for clipping displayed values.
#' @return ggplot point layer. Rank on the original y scale before clipping and
#'   transforming for display. Grid sampling restores the previous random state.
#' @noRd
stat_qq_points <- function(
    mapping = NULL,
    data = NULL,
    geom = "point",
    position = "identity",
    show.legend = NA,
    distribution = "unif",
    max_pts_to_plot = 500,
    ymin = -Inf,
    ymax = Inf,
    inherit.aes = TRUE,
    ...
) {
    ggplot2::layer(
        stat = StatQQPoints,
        data = data,
        mapping = mapping,
        geom = geom,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(
            distribution = distribution,
            max_pts_to_plot = max_pts_to_plot,
            ymin = ymin,
            ymax = ymax,
            ...
        )
    )
}


#' Transform a positive value to minus its logarithm for plotting
#' @param base Logarithm base, defaulting to e.
#' @return Scale transformation -log(x, base), with inverse base^(-x), logarithmic
#'   breaks and domain `[1e-100, Inf]`. The two local functions implement the forward
#'   and inverse calculations; values are not otherwise checked here.
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
