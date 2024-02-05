get_my_theme <- function(element_text_size = 11) {
  ggplot2::theme_bw() + ggplot2::theme(axis.line = ggplot2::element_line(color = "black"),
                                       panel.grid.major = ggplot2::element_blank(),
                                       panel.grid.minor = ggplot2::element_blank(),
                                       panel.border = ggplot2::element_blank(),
                                       panel.background = ggplot2::element_blank(),
                                       plot.title = ggplot2::element_text(hjust = 0.5, size = element_text_size))
}

######################
# 1. PLOT ASSIGN GRNAS
######################
#' Plot gRNA count distributions
#'
#' `plot_grna_count_distributions()` plots the empirical UMI count distribution of one or more gRNAs.
#'
#' @param sceptre_object any initialized `sceptre_object`. `plot_grna_count_distributions()` can be called at any point in the pipeline after `import_data()`.
#' @param n_grnas_to_plot (optional; default \code{4}) an integer specifying the number of randomly selected gRNAs to plot
#' @param grnas_to_plot (optional; default `NULL`) a character vector giving the names of one or more specific gRNAs to plot. If \code{NULL}, then `n_grnas_to_plot` random gRNAs are plotted.
#' @param threshold (optional; default `NULL`) an integer representing a gRNA count cut-off; if provided, the bins of length 1 will go up to and include this value, after which the exponentially growing bins begin. A vertical line is also drawn at this value. If \code{NULL}, then 10 is the largest gRNA count with its own bin. Non-integer values will be rounded.
#'
#' @details
#' The x-axis is a piecewise linear-log scale, with bins of size 1 going from gRNA counts of 0 up to max(10, \code{threshold}), and then the bin widths grow exponentially in size. The number under each bar indicates the first value that is counted for that bar, and that bar includes all integers from that label up until the integer immediately preceding the label of the next bar on the right. For example, if one bar has a label of "23" and the next bar on the right has a label of "26" then the bar with the label of "23" counts values of 23, 24, and 25 in the data.
#'
#' @return a single \code{ggplot2} plot
#' @export
#' @examples
#' # A full example can be found at ?sceptre
plot_grna_count_distributions <- function(sceptre_object, n_grnas_to_plot = 4L, grnas_to_plot = NULL, threshold = NULL) {
  grna_matrix <- get_grna_matrix(sceptre_object)
  # rounding just in case the user provides a non-integer one
  if (!is.null(threshold)) threshold <- round(threshold)
  if (is.null(grnas_to_plot)) {
    grnas_to_plot <- sample(x = rownames(grna_matrix), size = min(nrow(grna_matrix), n_grnas_to_plot), replace = FALSE)
  } else {
    if (!(all(grnas_to_plot %in% rownames(grna_matrix)))) stop("gRNA IDs must be a subset of the rownames of the gRNA matrix.")
  }
  grna_matrix <- set_matrix_accessibility(grna_matrix, make_row_accessible = TRUE)
  grna_expressions <- lapply(X = grnas_to_plot, function(grna_id) {
    load_row(grna_matrix, grna_id)
  }) |> unlist()
  grna_ids_rep <- rep(factor(grnas_to_plot), each = ncol(grna_matrix))
  to_plot <- data.frame(grna_id = grna_ids_rep, grna_expressions = grna_expressions) |>
    dplyr::filter(grna_expressions < 10000)

  # this function takes a vector of grna expressions and returns a data.frame
  # which gets passed to `cut` for binning that vector. In the returned data.frame,
  # one column (`bin_upper_bounds`) contains the upper end point of the resulting bin,
  # and the other column (`bin_labels`) has the name that that bin will get.
  grna_expressions_to_binned_factor <- function(gnra_expressions) {
    max_expression_count <- max(gnra_expressions)
    max_single_bin <- max(10, threshold - 1) # 10 is just a nice convenient default
    bin_upper_bounds <- 0:max_single_bin
    bin_labels <-  as.character(bin_upper_bounds)
    if (max_expression_count > 10) { # now we need to add exp growing bins, w/ more complex labels
      # this next section relies on the fact that the upper bounds are going to be at locations
      # max_single_bin + 2, max_single_bin + 2 + 2^2, max_single_bin + 2 + 2^2 + 2^3, ...
      # and 2 + 2^2 + 2^3 + ... + 2^n = 2(2^n-1), so `num_exp_bins` comes from finding the
      # smallest n such that this biggest bin width is above `max_expression_count`
      num_exp_bins <- log2((max_expression_count - max_single_bin) / 2 + 1) |> ceiling()
      bin_upper_bounds <- c(bin_upper_bounds, max_single_bin + 2 * (2^(1:num_exp_bins) - 1))
      bin_labels <- c(bin_labels, as.character(max_single_bin + (2^((1:num_exp_bins) - 1) - 1)))
    }
    return(data.frame(bin_upper_bounds = bin_upper_bounds, bin_labels = bin_labels))
  }

  # creating a list of bar plots for each grna_id, so that we can avoid dropping
  # any bins inside the range of the expressions for each grna_id, but we don't
  # keep empty ones in the tail.
  plot_list <- lapply(grnas_to_plot, function(curr_grna_id) {
    curr_df <- dplyr::filter(to_plot, as.character(grna_id) == curr_grna_id) |>
      dplyr::mutate(grna_id  = droplevels(grna_id))
    bin_info <- grna_expressions_to_binned_factor(curr_df$grna_expressions)
    p <- curr_df |>
      dplyr::mutate(
        grna_expressions_bin = cut(grna_expressions, breaks = c(-Inf, bin_info$bin_upper_bounds),
                                   labels = bin_info$bin_labels)
      ) |>
      dplyr::group_by(grna_id, grna_expressions_bin) |>
      dplyr::summarize(bin_counts = dplyr::n(), .groups = "drop_last") |>
      ggplot2::ggplot(mapping = ggplot2::aes(x = grna_expressions_bin , y = bin_counts)) +
      ggplot2::geom_bar(stat = "identity", fill = "grey90", col = "midnightblue") +
      ggplot2::scale_y_continuous(
        trans = scales::pseudo_log_trans(base = 10, sigma = 0.5),
        breaks = c(0, 10, 100, 1000, 100000), expand = c(0, NA)
      ) +
      get_my_theme() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
                     axis.title = ggplot2::element_blank()
      ) +
      ggplot2::ggtitle(curr_grna_id) +
      ggplot2::scale_x_discrete(drop = FALSE)
    if (!is.null(threshold)) {
      # adding +.5 so it is after the bin rather than through the middle
      p <- p + ggplot2::geom_vline(xintercept = threshold + .5, color = "mediumseagreen", linetype = "dashed")
    }
    return(p)
  })

  # these are the dimensions of our cowplot grid
  n_row <-  floor(sqrt(length(grnas_to_plot)))
  n_col <- ceiling(length(grnas_to_plot) / n_row)

  # for plots on the left-most edge, add a y-axis label
  for(i in 0:(n_row - 1))  {
    plot_list[[1 + i * n_col]] <- plot_list[[1 + i * n_col]] +
      ggplot2::ylab("Cell count (log scale)") +
      ggplot2::theme(axis.title.y = ggplot2::element_text())
  }
  # for plots on the bottom row, add a x-axis label
  for(j in (1 + (n_row - 1) * n_col):length(grnas_to_plot))  {
    plot_list[[j]] <- plot_list[[j]] + ggplot2::xlab("gRNA count") +
      ggplot2::theme(axis.title.x = ggplot2::element_text())
  }
  p <- do.call(
    what = cowplot::plot_grid,
    args = c(plot_list,  nrow = n_row)
  )
  return(p)
}

#' Plot assign gRNAs
#'
#' `plot_assign_grnas()` plots the outcome of the gRNA-to-cell assignment step. See \href{https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_assign_grnas}{Section 3 of the introductory chapter in the manual} for guidance on interpreting this visual.
#'
#' @param sceptre_object a \code{sceptre_object} that has had `assign_grnas()` called on it
#' @param n_grnas_to_plot (optional; default \code{3}) the number of gRNAs to display in the plots of gRNA count versus cell assignment
#' @param grnas_to_plot (optional; default \code{NULL}) a character vector giving the names of specific gRNAs to plot; if \code{NULL}, then the gRNAs are chosen at random.
#' @param point_size (optional; default \code{0.9}) the size of the individual points in the plot
#' @param transparency (optional; default \code{0.8}) the transparency of the individual points in the plot
#' @param return_indiv_plots (optional; default \code{FALSE}) if \code{FALSE}, then a list of \code{ggplot} objects is returned; if \code{TRUE} then a single \code{cowplot} object is returned.
#'
#' @return a single \code{cowplot} object containing the combined panels (if \code{return_indiv_plots} is set to \code{TRUE}) or a list of the individual panels (if \code{return_indiv_plots} is set to \code{FALSE})
#' @export
#' @examples
#' # A full example can be found at ?sceptre;
#' # `plot_assign_grnas()` is dispatched when
#' # `plot()` is called on the `sceptre_object`
#' # in step 3 (the gRNA assignment step).
plot_assign_grnas <- function(sceptre_object, n_grnas_to_plot = 3L, grnas_to_plot = NULL, point_size = 0.9, transparency = 0.8, return_indiv_plots = FALSE) {
  n_points_to_plot_per_umi <- 1000
  n_grnas_to_plot_panel_b <- 1000
  if (!sceptre_object@functs_called["assign_grnas"]) {
    stop("This `sceptre_object` has not yet had `assign_grnas` called on it.")
  }
  init_assignments <- sceptre_object@initial_grna_assignment_list
  grna_matrix <- get_grna_matrix(sceptre_object) |> set_matrix_accessibility(make_row_accessible = TRUE)
  grna_ids <- names(init_assignments)
  # sample grnas to plot
  if (is.null(grnas_to_plot)) {
    grnas_to_plot <- sample(x = grna_ids, size = min(nrow(grna_matrix), n_grnas_to_plot), replace = FALSE)
  } else {
    if (!(all(grnas_to_plot %in% grna_ids))) stop("gRNA IDs must be a subset of the rownames of the gRNA matrix.")
  }
  to_plot_a <- lapply(X = grnas_to_plot, function(curr_grna_id) {
    assignment <- multiple_grnas <- logical(length = ncol(grna_matrix)) # logical vecs w/ one entry per cell
    assignment[init_assignments[[curr_grna_id]]] <- TRUE # for this grna, `assignment` indicates which cells got this grna initially
    multiple_grnas[sceptre_object@cells_w_multiple_grnas] <- TRUE  # indicates which cells have >1 grna
    g <- if (nrow(sceptre_object@grna_target_data_frame_with_vector) >= 1L) {
      sceptre_object@grna_target_data_frame_with_vector |>
        dplyr::filter(vector_id == curr_grna_id) |>
        dplyr::pull(grna_id) |>
        sapply(function(i) load_row(grna_matrix, i)) |>
        rowSums()
    } else {
      load_row(grna_matrix, curr_grna_id)
    }
    df <- data.frame(g = g,
                     assignment = ifelse(assignment, "pert", "unpert") |> factor(),
                     grna_id = curr_grna_id |> factor(),
                     multiple_grnas = multiple_grnas)
    # if assignment method maximum, remove cells containing multiple gRNAs
    if (sceptre_object@grna_assignment_method == "maximum") df <- df |> dplyr::filter(!multiple_grnas)
    return(df)
  }) |> data.table::rbindlist()

  # downsample the unperturbed cells
  to_plot_a <- to_plot_a |>
    dplyr::group_by(g) |>
    dplyr::sample_n(size = min(n_points_to_plot_per_umi, dplyr::n())) |>
    dplyr::ungroup()

  # plot a
  p_a <- ggplot2::ggplot(data = to_plot_a, mapping = ggplot2::aes(x = assignment, y = g, col = assignment)) +
    ggplot2::geom_jitter(alpha = transparency, size = point_size) +
    ggplot2::facet_wrap(nrow = 1, facets = grna_id ~ ., scales = "free_y") +
    ggplot2::scale_y_continuous(trans = "log1p", breaks = c(0, 1, 5, 15, 50, 200, 1000, 5000, 10000)) +
    ggplot2::xlab("Assignment") +
    ggplot2::ylab("gRNA count") +
    get_my_theme() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
    ggplot2::scale_color_manual(values = c("firebrick1", "darkorchid1"))

  # plot b
  n_cells_per_grna <- sapply(init_assignments, length)
  mean_c_cells_per_grna <- mean(n_cells_per_grna)
  to_plot_b <- data.frame(x = n_cells_per_grna,
                          y = names(n_cells_per_grna)) |>
    dplyr::arrange(n_cells_per_grna) |>
    dplyr::mutate(y = factor(y, labels = y, levels = y)) |>
    dplyr::sample_n(min(n_grnas_to_plot_panel_b, dplyr::n()))
  p_b <- ggplot2::ggplot(data = to_plot_b,
                         mapping = ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_bar(stat = "identity", width = 0.5, fill = "grey90", col = "darkblue") +
    get_my_theme() + ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::xlab("N cells") + ggplot2::ylab("gRNA") +
    ggplot2::ggtitle("N cells per gRNA") +
    ggplot2::geom_vline(xintercept = mean_c_cells_per_grna, col = "darkorchid1", lwd = 1.0) +
    ggplot2::annotate(geom = "text", label = paste0("Mean cells per\ngRNA = ", round(mean_c_cells_per_grna, 2)),
                      x = Inf, y = -Inf, vjust = -0.5, hjust = 1, col = "darkorchid1", size = 3.0)

  # plot c
  if (sceptre_object@grna_assignment_method == "maximum") {
    p_c <- ggplot2::ggplot() +
      ggplot2::theme_void()
  } else {
    n_grnas_per_cell <- sceptre_object@grnas_per_cell
    moi <- mean(n_grnas_per_cell)
    p_c <- ggplot2::ggplot(data = data.frame(x = n_grnas_per_cell),
                           mapping = ggplot2::aes(x = x)) +
      ggplot2::geom_histogram(binwidth = max(1, 0.02 * length(unique(n_grnas_per_cell))),
                              fill = "grey90", color = "darkblue") +
      ggplot2::scale_y_continuous(expand = c(0, 0), trans = "log1p", breaks = 10^(1:8)) +
      get_my_theme() + ggplot2::ylab("Frequency") +
      ggplot2::ggtitle("N gRNAs per cell") +
      ggplot2::xlab("N gRNAs") +
      ggplot2::geom_vline(xintercept = mean(n_grnas_per_cell), col = "darkorchid1", lwd = 1.0) +
      ggplot2::annotate(geom = "text", label = paste0("MOI = ", round(moi, 2)),
                        x = Inf, y = Inf, vjust = 2.0, hjust = 1, col = "darkorchid1", size = 3.0)
  }
  if (return_indiv_plots) {
    out <- list(p_a, p_b, p_c)
  } else {
    bottom <- cowplot::plot_grid(p_b, p_c, nrow = 1)
    out <- cowplot::plot_grid(p_a, bottom, ncol = 1)
  }
  return(out)
}

###########################
# 3. PLOT CALIBRATION CHECK
###########################
#' Plot run calibration check
#'
#' \code{plot_run_calibration_check()} creates a visualization of the outcome of the calibration check. See \href{https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_calibration_check}{Section 5 of the introductory chapter in the manual} for guidance on interpreting this visual.
#'
#' @param sceptre_object a \code{sceptre_object} that has had \code{run_calibration_check} called on it
#' @param point_size (optional; default \code{0.55}) the size of the individual points in the plot
#' @param transparency (optional; default \code{0.8}) the transparency of the individual points in the plot
#' @param return_indiv_plots (optional; default \code{FALSE}) if \code{FALSE} then a list of \code{ggplot} is returned; if \code{TRUE} then a single \code{cowplot} object is returned.
#'
#' @return a single \code{cowplot} object containing the combined panels (if \code{return_indiv_plots} is set to \code{TRUE}) or a list of the individual panels (if \code{return_indiv_plots} is set to \code{FALSE})
#' @export
#'
#' @examples
#' # A full example can be found at ?sceptre;
#' # `plot_run_calibration_check()` is dispatched when
#' # `plot()` is called on the `sceptre_object`
#' # in step 5 (the run calibration check step).
plot_run_calibration_check <- function(sceptre_object, point_size = 0.55, transparency = 0.8, return_indiv_plots = FALSE) {
  if (!sceptre_object@functs_called["run_calibration_check"]) {
    stop("This `sceptre_object` has not yet had `run_calibration_check` called on it.")
  }
  calibration_result <- sceptre_object@calibration_result
  my_theme <- get_my_theme()

  # compute n rejections
  n_rejections <- sum(calibration_result$significant)
  str1 <- paste0("Number of\nfalse discoveries\n(at alpha ", signif(sceptre_object@multiple_testing_alpha, 1), "): ", n_rejections)
  str2 <- paste0("\n\nMean log-2-fold\nchange: ", signif(mean(calibration_result$log_2_fold_change), 2))
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
                         mapping = ggplot2::aes(x = log_2_fold_change)) +
    ggplot2::geom_histogram(binwidth = if (nrow(calibration_result) > 10000) 0.02 else 0.05,
                            fill = "grey90", col = "black", boundary = 0) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.0, .01))) +
    ggplot2::ggtitle("Log fold changes") +
    ggplot2::xlab("Estimated log-2 fold change") + ggplot2::ylab("Density") +
    ggplot2::geom_vline(xintercept = 0, col = "firebrick2", linewidth = 1) +
    my_theme

  p_d <- ggplot2::ggplot() +
    ggplot2::annotate(geom = "text", label = str, x = 1.3, y = 1.2) +
    ggplot2::theme_void() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", color = "white")) +
    ggplot2::xlim(c(0, 2)) +
    ggplot2::ylim(c(0, 2))

  if (return_indiv_plots) {
    out <- list(p_a, p_b, p_c, p_d)
  } else {
    out <- cowplot::plot_grid(p_a, p_b, p_c, p_d, nrow = 2, rel_heights = c(0.55, 0.45))
  }
  return(out)
}

############################
# 4. PLOT DISCOVERY ANALYSIS
############################
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
      ggplot2::ggtitle("QQ plot (tail)") +
      (if (!is.na(p_thresh)) ggplot2::geom_hline(yintercept = p_thresh, linetype = "dashed") else NULL)
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
    ggplot2::scale_y_continuous(trans = revlog_trans(10), expand = c(0.02, 0)) +
    get_my_theme() + ggplot2::xlab("Log fold change") + ggplot2::ylab("P-value") +
    (if (!is.na(p_thresh)) ggplot2::geom_hline(yintercept = p_thresh, linetype = "dashed") else NULL) +
    ggplot2::theme(legend.position = "none") + ggplot2::scale_color_manual(values = c("dodgerblue", "blueviolet")) +
    ggplot2::ggtitle("Discovery volcano plot")
  return(out)
}


#' Plot run discovery analysis
#'
#' `plot_run_discovery_analysis()` creates a visualization of the outcome of the discovery analysis. See \href{https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_run_discovery_analysis}{Section 7 of the introductory chapter in the manual} for guidance on interpreting this visual.
#'
#' @param sceptre_object a \code{sceptre_object} that has had \code{run_discovery_analysis} called on it
#' @param x_limits (optional; default \code{c(-1.5, 1.5)}) a numeric vector of length 2 giving the lower and upper limits of the x-axis (corresponding to log-2 fold change) for the "Discovery volcano plot" panel
#' @param point_size (optional; default \code{0.55}) the size of the individual points in the plot
#' @param transparency (optional; default \code{0.8}) the transparency of the individual points in the plot
#' @param return_indiv_plots (optional; default \code{FALSE}) if \code{FALSE} then a list of \code{ggplot} is returned; if \code{TRUE} then a single \code{cowplot} object is returned.
#' @return  a single \code{cowplot} object containing the combined panels (if \code{return_indiv_plots} is set to \code{TRUE}) or a list of the individual panels (if \code{return_indiv_plots} is set to \code{FALSE})
#'
#' @export
#' @examples
#' # A full example can be found at ?sceptre;
#' # `plot_run_discovery_analysis()` is dispatched when
#' # `plot()` is called on the `sceptre_object`
#' # in step 7 (the run discovery analysis step).
plot_run_discovery_analysis <- function(sceptre_object, x_limits = c(-1.5, 1.5), point_size = 0.55, transparency = 0.8, return_indiv_plots = FALSE) {
  if (!sceptre_object@functs_called["run_discovery_analysis"]) {
    stop("This `sceptre_object` has not yet had `run_discovery_analysis` called on it.")
  }
  # first, compute the rejection set
  discovery_result <- sceptre_object@discovery_result |> stats::na.omit()
  calibration_result <- sceptre_object@calibration_result
  discovery_set <- discovery_result |> dplyr::filter(significant)
  p_thresh <- if (nrow(discovery_set) >= 1L) max(discovery_set$p_value) else NA
  if (nrow(calibration_result) != nrow(discovery_result)) { # if the two sets of pairs do not coincide in number, set calibration_result to data.frame()
    calibration_result <- data.frame()
  }
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
  discovery_result_downsample <- downsample_result_data_frame(result_df = discovery_result)
  p3 <- make_volcano_plot(discovery_result = discovery_result_downsample,
                          p_thresh = p_thresh,
                          transparency = transparency,
                          point_size = point_size,
                          x_limits = x_limits)
  # create the final panel
  n_rejections <- nrow(discovery_set)
  n_pairs <- nrow(discovery_result)
  str <- paste0("Number of discovery pairs\ncalled as significant\n(at alpha ", signif(sceptre_object@multiple_testing_alpha, 1), "): ", n_rejections, " of ", n_pairs)
  p4 <- ggplot2::ggplot() +
    ggplot2::annotate(geom = "text", label = str, x = 1.1, y = 1.2) +
    ggplot2::theme_void() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", color = "white")) +
    ggplot2::xlim(c(0, 2)) +
    ggplot2::ylim(c(0, 2))

  # combine the plots
  if (return_indiv_plots) {
    p_out <- list(p1, p2, p3, p4)
  } else {
    p_out <- cowplot::plot_grid(p1, p2, p3, p4, labels = c("a", "b", "c"), rel_heights = c(0.55, 0.45), nrow = 2)
  }
  return(p_out)
}

############
# 5. PLOT QC
############
#' Plot covariates
#'
#' `plot_covariates()` creates histograms of several cell-specific covariates: `response_n_nonzero`, `response_n_umis`, and (if applicable) `response_p_mito`.
#'
#' To help guide the selection of QC thresholds, `plot_covariates()` plots candidate QC thresholds as vertical lines on the histograms. If `run_qc()` has been called on the `sceptre_object`, then the QC thresholds plotted on the histograms are those contained within the `sceptre_object`.
#'
#' @param sceptre_object any initialized `sceptre_object`. `plot_covariates()` can be called at any point in the pipeline after \code{import_data()}.
#' @param response_n_umis_range (optional; default \code{c(0.01, 0.99)}) a length-2 vector of quantiles indicating the location at which to draw vertical lines on the `response_n_umis` histogram
#' @param response_n_nonzero_range (optional; default \code{c(0.01, 0.99)}) a length-2 vector of quantiles indicating the location at which to draw vertical lines on the `response_n_nonzero` histogram
#' @param p_mito_threshold (optional; default \code{0.2}) a single numeric value in the interval \[0,1\] specifying the location at which to draw a vertical line on the `response_p_mito` histogram. Note that `p_mito_threshold` is an absolute number rather than a percentile.
#' @param return_indiv_plots (optional; default \code{FALSE}) if \code{FALSE}, then a list of \code{ggplot} objects is returned; if \code{TRUE} then a single \code{cowplot} object is returned.
#'
#' @return a single \code{cowplot} object containing the combined panels (if \code{return_indiv_plots} is set to \code{TRUE}) or a list of the individual panels (if \code{return_indiv_plots} is set to \code{FALSE})
#' @export
#' @examples
#' # A full example can be found at ?sceptre
plot_covariates <- function(sceptre_object,
                            response_n_umis_range = c(0.01, 0.99),
                            response_n_nonzero_range = c(0.01, 0.99),
                            p_mito_threshold = 0.2,
                            return_indiv_plots = FALSE) {
  if (sceptre_object@functs_called[["run_qc"]]) {
    response_n_umis_range <- sceptre_object@cellwise_qc_thresholds$response_n_umis_range
    response_n_nonzero_range <- sceptre_object@cellwise_qc_thresholds$response_n_nonzero_range
    p_mito_threshold <- sceptre_object@cellwise_qc_thresholds$p_mito_threshold
  }
  covariate_data_frame <- sceptre_object@covariate_data_frame
  make_histogram <- function(v, curr_range, plot_tit, use_quantile) {
    cutoffs <- if (use_quantile) stats::quantile(v, probs = curr_range) else curr_range
    p1 <- ggplot2::ggplot(data = data.frame(x = v),
                          mapping = ggplot2::aes(x = x)) +
      ggplot2::geom_histogram(col = "darkblue", fill = "grey90", bins = 50) +
      get_my_theme() +
      ggplot2::scale_y_continuous(expand = c(0, NA)) +
      ggplot2::geom_vline(xintercept = cutoffs[1], col = "darkorchid1", lwd = 1.0) +
      (if (length(cutoffs) == 2) {
        ggplot2::geom_vline(xintercept = cutoffs[2], col = "darkorchid1", lwd = 1.0)
      } else NULL) +
      ggplot2::ggtitle(plot_tit) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank())
  }
  p1 <- make_histogram(covariate_data_frame$response_n_nonzero, response_n_nonzero_range,
                       "Response N nonzero", use_quantile = TRUE)
  p2 <- make_histogram(covariate_data_frame$response_n_umis, response_n_umis_range,
                       "Response N UMIs", use_quantile = TRUE)
  p_mito_present <- "response_p_mito" %in% colnames(covariate_data_frame)
  if (p_mito_present) {
    p3 <- make_histogram(covariate_data_frame$response_p_mito, p_mito_threshold,
                         plot_tit = "Percent mito", use_quantile = FALSE)
    p_out <- if (return_indiv_plots) list(p1, p2, p3) else cowplot::plot_grid(p1, p2, p3, NULL, ncol = 2)
  } else {
    p_out <- if (return_indiv_plots) list(p1, p2) else cowplot::plot_grid(p1, p2, ncol = 2)
  }
  return(p_out)
}


#' Plot run QC
#'
#' `plot_run_qc()` creates a visualization of the outcome of the QC step. See \href{https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_qc}{Section 4 of the introductory chapter in the manual} for guidance on interpreting this visual.
#'
#' @param sceptre_object a \code{sceptre_object} that has had \code{run_qc()} called on it
#' @param downsample_pairs (optional; default \code{10000}) the maximum number of points to plot in the lower panel of the figure (i.e., the "pairwise QC" plot)
#' @param point_size (optional; default \code{0.55}) the size of the individual points in the plot
#' @param transparency (optional; default \code{0.8}) the transparency of the individual points in the plot
#' @param return_indiv_plots (optional; default \code{FALSE}) if \code{FALSE} then a list of \code{ggplot} is returned; if \code{TRUE} then a single \code{cowplot} object is returned.
#' @return  a single \code{cowplot} object containing the combined panels (if \code{return_indiv_plots} is set to \code{TRUE}) or a list of the individual panels (if \code{return_indiv_plots} is set to \code{FALSE})
#'
#' @export
#' @examples
#' # A full example can be found at ?sceptre;
#' # `plot_run_qc()` is dispatched when
#' # `plot()` is called on the `sceptre_object`
#' # in step 4 (the run qc step).
plot_run_qc <- function(sceptre_object, downsample_pairs = 10000L, point_size = 0.55, transparency = 0.8, return_indiv_plots = FALSE) {
  if (!sceptre_object@functs_called["run_qc"]) {
    stop("This `sceptre_object` has not yet had `run_qc` called on it.")
  }
  my_cols <- c("mediumseagreen", "indianred2")
  my_breaks <- c(0, 1, 3, 7, 50, 500, 5000, 50000)
  discovery_pairs <- sceptre_object@discovery_pairs_with_info |>
    dplyr::mutate(pass_qc = ifelse(pass_qc, "Pass", "Fail")) |>
    dplyr::mutate(pass_qc = factor(pass_qc, levels = c("Pass", "Fail")))
  if (nrow(discovery_pairs) >= downsample_pairs) discovery_pairs <- discovery_pairs |> dplyr::sample_n(downsample_pairs)

  # pairwise QC
  p_b <- ggplot2::ggplot(data = discovery_pairs,
                  mapping = ggplot2::aes(x = n_nonzero_trt, y = n_nonzero_cntrl, col = pass_qc)) +
    ggplot2::geom_point(alpha = transparency, size = point_size) + get_my_theme() +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 1),
                                breaks = my_breaks) +
    ggplot2::scale_x_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 1),
                                breaks = my_breaks) +
    ggplot2::geom_hline(yintercept = sceptre_object@n_nonzero_cntrl_thresh) +
    ggplot2::geom_vline(xintercept = sceptre_object@n_nonzero_trt_thresh) +
    ggplot2::xlab("N nonzero trt. cells") +
    ggplot2::ylab("N nonzero cntrl. cells") +
    ggplot2::theme(legend.position = "bot") +
    ggplot2::scale_color_manual(values = my_cols) +
    ggplot2::ggtitle("Pairwise QC") +
    ggplot2::theme(legend.position = c(0.85, 0.2),
                   legend.margin = ggplot2::margin(t = -0.5, unit = "cm"),
                   legend.title = ggplot2::element_blank()) +
    ggplot2::guides(color = ggplot2::guide_legend(
      keywidth = 0.0,
      keyheight = 0.15,
      default.unit = "inch",
      override.aes = list(size = 1.25)))

  # cellwise QC
  cell_removal_metrics <- sceptre_object@cell_removal_metrics
  n_orig_cells <- ncol(get_response_matrix(sceptre_object))
  cell_removal_metrics_frac <- cell_removal_metrics/n_orig_cells * 100
  df <- data.frame(fraction_cells_removed = cell_removal_metrics_frac,
                   Filter = c("N response UMIs", "N nonzero responses", "Percent mito", "multiple gRNAs", "User-specified", "Any filter"))
  if (!sceptre_object@low_moi) df <- df |> dplyr::filter(Filter != "multiple gRNAs")
  # make a barplot. remove x-axis text
  p_a <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = Filter, y = fraction_cells_removed)) +
    ggplot2::geom_bar(stat = "identity", fill = "grey90", col = "darkblue") + get_my_theme() +
    ggplot2::scale_y_continuous(expand = c(0, NA)) +
    ggplot2::ylab("Percent cells removed") +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 20)) +
    get_my_theme() +
    ggplot2::theme(legend.position = "none", axis.title.x = ggplot2::element_blank()) +
    ggplot2::ggtitle("Cellwise QC")

  # combine the plots
  if (return_indiv_plots) {
    p_out <- list(p_a, p_b)
  } else {
    p_out <- cowplot::plot_grid(p_a, p_b, ncol = 1)
  }
  return(p_out)
}

###############
# 6. PLOT POWER
###############
#' Plot run power check
#'
#' \code{plot_run_power_check()} creates a visualization of the outcome of the power check analysis. Positive (resp., negative) control p-values are plotted in the left (resp., right) column on a negative log scale.
#'
#' @param sceptre_object a \code{sceptre_object} that has had \code{run_power_check()} called on it
#' @param point_size (optional; default \code{1}) the size of the individual points in the plot
#' @param transparency (optional; default \code{0.8}) the transparency of the individual points in the plot
#' @param clip_to (optional; default \code{1e-20}) p-values smaller than this value are set to \code{clip_to} for better visualization. If \code{clip_to=0} is used then no clipping is done.
#'
#' @return a single \code{ggplot2} plot
#' @export
#' @examples
#' # A full example can be found at ?sceptre;
#' # `plot_run_power_check()` is dispatched when
#' # `plot()` is called on the `sceptre_object`
#' # in step 6 (the run power check step).
plot_run_power_check <- function(sceptre_object, point_size = 1, transparency = 0.8, clip_to = 1e-20) {
  if (!sceptre_object@functs_called["run_power_check"]) {
    stop("This `sceptre_object` has not yet had `run_power_check` called on it.")
  }
  my_theme <- get_my_theme()
  set.seed(3)
  my_cols <- c("mediumseagreen", "firebrick1")

  pos_ctrl_pvals <- sceptre_object@power_result$p_value |> stats::na.omit()
  neg_ctrl_pval_sub <- downsample_result_data_frame(
      result_df = sceptre_object@calibration_result
    ) |>
    dplyr::pull(p_value)
  group_names <- c("Positive Control", "Negative Control")
  df <- data.frame(
    lab = rep(
      group_names,
      c(length(pos_ctrl_pvals), length(neg_ctrl_pval_sub))
    ) |>
      factor(levels = group_names),
    p_values = c(pos_ctrl_pvals, neg_ctrl_pval_sub) |>
      pmax(clip_to)
  )

  p <- ggplot2::ggplot(
    data = df,
    mapping = ggplot2::aes(x = lab, y = p_values, color = lab)
  )
  my_breaks <- 10^(seq(-2, -40, by = -4))
  p <- p +
    ggplot2::geom_jitter(width = .25, height = 0, size = point_size, alpha = transparency) +
    ggplot2::scale_y_continuous(trans = revlog_trans(base = 10), expand = c(0.01, 0), breaks = my_breaks) +
    ggplot2::scale_color_manual(values = my_cols, guide = "none") +
    ggplot2::labs(
      x = "Pair type",
      y = "p-value",
      title = "Positive and negative control p-values") +
    my_theme

  return(p)
}


downsample_result_data_frame <- function(result_df, downsample_pairs = 1000) {
  result_df |>
    dplyr::mutate(bin = cut(p_value, breaks = c(10^(-seq(0, 20)), 0))) |>
    dplyr::group_by(bin) |>
    dplyr::sample_n(size = min(dplyr::n(), downsample_pairs))
}
