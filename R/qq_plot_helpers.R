# workhorse function for \code{stat_qq_band}
StatQQBand <- ggplot2::ggproto(
    "StatQQBand",
    ggplot2::Stat,
    required_aes = c("y"), # the sample will come through the y aesthetic
    #' Select one equal-sized QQ group per panel before band computation.
    #' @param data Transformed y, PANEL and ggplot-assigned group IDs.
    #' @param params Existing ggplot stat parameters.
    #' @return Input rows for the minimum group separately in each panel, retaining
    #'   row order; errors on unequal group sizes within a panel. No sampling/RNG.
    #' @noRd
    setup_data = function(data, params) {
        if ("group" %in% names(data)) {
            # Require equal group sizes within each panel before selecting a band.
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
    #' Calculate the existing pointwise order-statistic QQ band for one group.
    #' @param data,scales Transformed y and panel scales; y is inverted first.
    #' @param distribution,max_pts_to_plot,ci_level Quantile family, grid resolution
    #'   and pointwise confidence. Rank ties use first occurrence.
    #' @return Transformed x/y/ymin/ymax rows, one per occupied x/y grid cell under
    #'   seed 4. Restores caller RNG; no global point cap or new band method.
    #' @noRd
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

#' Add confidence band for a QQ plot
#'
#' @inheritParams ggplot2::stat_identity
#' @param distribution The distribution with respect to which to make QQ plot
#' (e.g. \code{"unif"}, \code{"norm"}, etc.); default \code{"unif"}
#' @param max_pts_to_plot Number of intervals per axis for deterministic x/y-grid
#' downsampling; at most one point is kept per occupied grid cell, not a global cap.
#' @param ci_level The pointwise confidence level for the QQ plot; default 0.95
#' @return A ggplot layer using the existing per-panel pointwise beta bands,
#' inverse/forward scales, first-occurrence ranks and restored seed-4 RNG.
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
    #' Calculate existing QQ points for one group in its panel scales.
    #' @param data,scales Transformed y and panel scales; ranks use inverse-scale y.
    #' @param distribution,max_pts_to_plot Quantile family and occupied-grid resolution.
    #' @param ymin,ymax Inverse-scale display clipping bounds, applied after ranking.
    #' @return Transformed x/y rows sampled per occupied grid cell with seed 4;
    #'   first-occurrence tie ranks and caller RNG are preserved.
    #' @noRd
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

#' Add points for a QQ plot
#'
#' @inheritParams ggplot2::stat_identity
#' @param distribution The distribution with respect to which to make QQ plot
#' (e.g. \code{"unif"}, \code{"norm"}, etc.); default \code{"unif"}
#' @param max_pts_to_plot Number of intervals per axis for deterministic x/y-grid
#' downsampling; at most one point is kept per occupied grid cell, not a global cap.
#' @param ymin Samples with value less than \code{ymin} will be truncated at
#' \code{ymin} for plotting.
#' @param ymax Samples with value greater than \code{ymax} will be truncated at
#' \code{ymax} for plotting.
#' @return A ggplot point layer preserving inverse-scale ranks, display clipping,
#' occupied-grid sampling and caller RNG.
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


#' Reverse logarithm transformation
#'
#' @param base The base with respect to which to take the logarithm (defaults to
#' e).
#'
#' @return A scale transformation with -log(x, base), inverse base^(-x), log
#' breaks and the existing positive domain. The two local transform closures share
#' this contract; no validation, mutation, sampling or files are involved.
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
