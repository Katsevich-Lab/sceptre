get_my_theme <- function(element_text_size = 11) {
  ggplot2::theme_bw() + ggplot2::theme(
    axis.line = ggplot2::element_line(color = "black"),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(hjust = 0.5, size = element_text_size)
  )
}

######################
# 1. PLOT ASSIGN GRNAS
######################
#' Plot gRNA count distributions
#'
#' `plot_grna_count_distributions()` plots the empirical UMI count distribution of one or more gRNAs. `plot_grna_count_distributions()` can be called on a `sceptre_object` at any point in the pipeline after `import_data()`.
#'
#' @param sceptre_object a `sceptre_object`
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
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' ) |> plot_grna_count_distributions()
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
    bin_labels <- as.character(bin_upper_bounds)
    if (max_expression_count > 10) { # now we need to add exp growing bins, w/ more complex labels
      # this next section relies on the fact that the upper bounds are going to be at locations
      # max_single_bin + 2, max_single_bin + 2 + 2^2, max_single_bin + 2 + 2^2 + 2^3, ...
      # and 2 + 2^2 + 2^3 + ... + 2^n = 2(2^n-1), so `num_exp_bins` comes from finding the
      # smallest n such that this biggest bin width is above `max_expression_count`
      num_exp_bins <- log2((max_expression_count - max_single_bin) / 2 + 1) |> ceiling()
      bin_upper_bounds <- c(bin_upper_bounds, max_single_bin + 2 * (2^seq_len(num_exp_bins) - 1))
      bin_labels <- c(bin_labels, as.character(max_single_bin + (2^((seq_len(num_exp_bins)) - 1) - 1)))
    }
    return(data.frame(bin_upper_bounds = bin_upper_bounds, bin_labels = bin_labels))
  }

  # creating a list of bar plots for each grna_id, so that we can avoid dropping
  # any bins inside the range of the expressions for each grna_id, but we don't
  # keep empty ones in the tail.
  plot_list <- lapply(grnas_to_plot, function(curr_grna_id) {
    curr_df <- dplyr::filter(to_plot, as.character(grna_id) == curr_grna_id) |>
      dplyr::mutate(grna_id = droplevels(grna_id))
    bin_info <- grna_expressions_to_binned_factor(curr_df$grna_expressions)
    p <- curr_df |>
      dplyr::mutate(
        grna_expressions_bin = cut(grna_expressions,
          breaks = c(-Inf, bin_info$bin_upper_bounds),
          labels = bin_info$bin_labels
        )
      ) |>
      dplyr::group_by(grna_id, grna_expressions_bin) |>
      dplyr::summarize(bin_counts = dplyr::n(), .groups = "drop_last") |>
      ggplot2::ggplot(mapping = ggplot2::aes(x = grna_expressions_bin, y = bin_counts)) +
      ggplot2::geom_bar(stat = "identity", fill = "grey90", col = "midnightblue") +
      ggplot2::scale_y_continuous(
        trans = scales::pseudo_log_trans(base = 10, sigma = 0.5),
        breaks = c(0, 10, 100, 1000, 100000), expand = c(0, NA)
      ) +
      get_my_theme() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
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
  n_row <- floor(sqrt(length(grnas_to_plot)))
  n_col <- ceiling(length(grnas_to_plot) / n_row)

  # for plots on the left-most edge, add a y-axis label
  for (i in 0:(n_row - 1)) {
    plot_list[[1 + i * n_col]] <- plot_list[[1 + i * n_col]] +
      ggplot2::ylab("Cell count (log scale)") +
      ggplot2::theme(axis.title.y = ggplot2::element_text())
  }
  # for plots on the bottom row, add a x-axis label
  for (j in (1 + (n_row - 1) * n_col):length(grnas_to_plot)) {
    plot_list[[j]] <- plot_list[[j]] + ggplot2::xlab("gRNA count") +
      ggplot2::theme(axis.title.x = ggplot2::element_text())
  }
  p <- do.call(
    what = cowplot::plot_grid,
    args = c(plot_list, nrow = n_row)
  )
  return(p)
}

#' Plot assign gRNAs
#'
#' `plot_assign_grnas()` plots the outcome of the gRNA-to-cell assignment step. The top panel plots the gRNA-to-cell assignments of `n_grnas_to_plot` (default 3) randomly selected gRNAs. In each plot the points represent cells; the vertical axis indicates the UMI count of the gRNA in a given cell, and the horizontal axis indicates whether the cell has been classified as “perturbed” (i.e., it contains the gRNA) or unperturbed (i.e., it does not contain the gRNA). Perturbed (resp., unperturbed) cells are shown in the left (resp., right) column. The bottom left panel is a barplot of the number of cells to which each gRNA has been mapped. Finally, the bottom right panel is a histogram of the number of gRNAs contained in each cell. The mean number of gRNAs per cell --- i.e., the MOI --- is displayed in purple text.
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
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' ) |>
#'   set_analysis_parameters() |>
#'   assign_grnas(method = "thresholding") |>
#'   plot_assign_grnas()
plot_assign_grnas <- function(sceptre_object, n_grnas_to_plot = 3L, grnas_to_plot = NULL, point_size = 0.9, transparency = 0.8, return_indiv_plots = FALSE) {
  n_points_to_plot_per_umi <- 1000
  n_grnas_to_plot_panel_b <- 1000
  if (!sceptre_object@functs_called["assign_grnas"]) {
    stop("This `sceptre_object` has not yet had `assign_grnas` called on it.")
  }
  init_assignments <- sceptre_object@initial_grna_assignment_list
  grna_matrix <- get_grna_matrix(sceptre_object) |> set_matrix_accessibility(make_row_accessible = TRUE)
  grna_ids <- names(init_assignments)
  assigned <- vapply(init_assignments, length, FUN.VALUE = integer(1)) >= 1
  grna_ids <- grna_ids[assigned]
  # sample grnas to plot
  if (is.null(grnas_to_plot)) {
    grnas_to_plot <- sample(x = grna_ids, size = min(length(grna_ids), n_grnas_to_plot), replace = FALSE)
  } else {
    if (!(all(grnas_to_plot %in% grna_ids))) stop("gRNA IDs must be a subset of the rownames of the gRNA matrix.")
  }
  to_plot_a <- lapply(X = grnas_to_plot, function(curr_grna_id) {
    assignment <- cells_w_zero_or_twoplus_grnas <- logical(length = ncol(grna_matrix)) # logical vecs w/ one entry per cell
    assignment[init_assignments[[curr_grna_id]]] <- TRUE # for this grna, `assignment` indicates which cells got this grna initially
    cells_w_zero_or_twoplus_grnas[sceptre_object@cells_w_zero_or_twoplus_grnas] <- TRUE # indicates which cells have 0/2+ grnas
    g <- load_row(grna_matrix, curr_grna_id)
    df <- data.frame(
      g = g,
      assignment = ifelse(assignment, "pert", "unpert") |> factor(),
      grna_id = curr_grna_id |> factor(),
      cells_w_zero_or_twoplus_grnas = cells_w_zero_or_twoplus_grnas
    )
    # if assignment method maximum, remove cells with 0/2+ grnas
    if (sceptre_object@grna_assignment_method == "maximum") df <- df |> dplyr::filter(!cells_w_zero_or_twoplus_grnas)
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
  n_cells_per_grna <- vapply(init_assignments, length, FUN.VALUE = integer(1))
  mean_c_cells_per_grna <- mean(n_cells_per_grna)
  to_plot_b <- data.frame(
    x = n_cells_per_grna,
    y = names(n_cells_per_grna)
  ) |>
    dplyr::arrange(n_cells_per_grna) |>
    dplyr::mutate(y = factor(y, labels = y, levels = y)) |>
    dplyr::sample_n(min(n_grnas_to_plot_panel_b, dplyr::n()))
  p_b <- ggplot2::ggplot(
    data = to_plot_b,
    mapping = ggplot2::aes(x = x, y = y)
  ) +
    ggplot2::geom_bar(stat = "identity", width = 0.5, fill = "grey90", col = "darkblue") +
    get_my_theme() +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    ) +
    ggplot2::xlab("N cells") +
    ggplot2::ylab("gRNA") +
    ggplot2::ggtitle("N cells per gRNA") +
    ggplot2::geom_vline(xintercept = mean_c_cells_per_grna, col = "darkorchid1", lwd = 1.0) +
    ggplot2::annotate(
      geom = "text", label = paste0("Mean cells per\ngRNA = ", round(mean_c_cells_per_grna, 2)),
      x = Inf, y = -Inf, vjust = -0.5, hjust = 1, col = "darkorchid1", size = 3.0
    )

  # plot c
  if (sceptre_object@grna_assignment_method == "maximum") {
    p_c <- ggplot2::ggplot() +
      ggplot2::theme_void()
  } else {
    n_grnas_per_cell <- sceptre_object@grnas_per_cell
    moi <- mean(n_grnas_per_cell)
    p_c <- ggplot2::ggplot(
      data = data.frame(x = n_grnas_per_cell),
      mapping = ggplot2::aes(x = x)
    ) +
      ggplot2::geom_histogram(
        binwidth = max(1, 0.02 * length(unique(n_grnas_per_cell))),
        fill = "grey90", color = "darkblue"
      ) +
      ggplot2::scale_y_continuous(expand = c(0, 0), trans = "log1p", breaks = 10^(seq(1, 8))) +
      get_my_theme() +
      ggplot2::ylab("Frequency") +
      ggplot2::ggtitle("N gRNAs per cell") +
      ggplot2::xlab("N gRNAs") +
      ggplot2::geom_vline(xintercept = mean(n_grnas_per_cell), col = "darkorchid1", lwd = 1.0) +
      ggplot2::annotate(
        geom = "text", label = paste0("MOI = ", round(moi, 2)),
        x = Inf, y = Inf, vjust = 2.0, hjust = 1, col = "darkorchid1", size = 3.0
      )
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
#' \code{plot_run_calibration_check()} creates a visualization of the outcome of the calibration check. The visualization consists of four panels, which we describe below.
#' -   The upper left panel is a QQ plot of the p-values plotted on an untransformed scale. The p-values ideally should lie along the diagonal line, indicating uniformity of the p-values in the *bulk* of the distribution.
#' -   The upper right panel is a QQ plot of the p-values plotted on a negative log-10 transformed scale. The p-values ideally should lie along the diagonal line (with the majority of the p-values falling within the gray confidence band), indicating uniformity of the p-values in the *tail* of the distribution.
#' -   The lower left panel is a histogram of the estimated log-2 fold changes. The histogram ideally should be roughly symmetric and centered around zero.
#' -   Finally, the bottom right panel is a text box displaying (i) the number of false discoveries that `sceptre` has made on the negative control data and (ii) the mean estimated log-fold change.
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
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' # import data
#' sceptre_object <- import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' )
#' sceptre_object |>
#'   set_analysis_parameters(
#'     side = "left",
#'     resampling_mechanism = "permutations"
#'   ) |>
#'   assign_grnas(method = "thresholding") |>
#'   run_qc() |>
#'   run_calibration_check(
#'     parallel = TRUE,
#'     n_processors = 2,
#'     n_calibration_pairs = 500,
#'     calibration_group_size = 2,
#'   ) |>
#'   plot_run_calibration_check()
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

  p_a <- ggplot2::ggplot(
    data = calibration_result,
    mapping = ggplot2::aes(y = p_value)
  ) +
    stat_qq_points(
      ymin = 1e-8, size = point_size,
      col = "firebrick2",
      alpha = transparency
    ) +
    stat_qq_band() +
    ggplot2::scale_x_reverse() +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(x = "Expected null p-value", y = "Observed p-value") +
    ggplot2::geom_abline(col = "black") +
    ggplot2::ggtitle("QQ plot (bulk)") +
    my_theme

  p_b <- ggplot2::ggplot(data = calibration_result, mapping = ggplot2::aes(y = p_value)) +
    stat_qq_points(
      ymin = 1e-8, size = point_size,
      col = "firebrick2", alpha = transparency
    ) +
    stat_qq_band() +
    ggplot2::scale_x_continuous(trans = revlog_trans(10)) +
    ggplot2::scale_y_continuous(trans = revlog_trans(10)) +
    ggplot2::labs(x = "Expected null p-value", y = "Observed p-value") +
    ggplot2::geom_abline(col = "black") +
    ggplot2::ggtitle("QQ plot (tail)") +
    my_theme

  p_c <- ggplot2::ggplot(
    data = calibration_result |> dplyr::filter(abs(log_2_fold_change) < 0.6),
    mapping = ggplot2::aes(x = log_2_fold_change)
  ) +
    ggplot2::geom_histogram(
      binwidth = if (nrow(calibration_result) > 10000) 0.02 else 0.05,
      fill = "grey90", col = "black", boundary = 0
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.0, .01))) +
    ggplot2::ggtitle("Log fold changes") +
    ggplot2::xlab("Estimated log-2 fold change") +
    ggplot2::ylab("Density") +
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
  lab <- c(
    rep(factor("Negative control"), nrow(calibration_result)),
    rep(factor("Discovery"), nrow(discovery_result))
  ) |>
    factor(levels = c("Discovery", "Negative control"))
  df <- data.frame(p_value = c(calibration_result$p_value, discovery_result$p_value), lab = lab)

  p_out <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(y = p_value, col = lab)) +
    stat_qq_points(ymin = 1e-8, size = point_size, alpha = transparency) +
    stat_qq_band() +
    ggplot2::labs(x = "Expected null p-value", y = "Observed p-value") +
    ggplot2::geom_abline(col = "black") +
    get_my_theme() +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      legend.position = if (include_legend) "bottom" else "none",
      axis.title.y = if (include_y_axis_text) ggplot2::element_text() else ggplot2::element_blank()
    ) +
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
      ggplot2::theme(
        legend.position.inside = c(0.7, 0.1),
        legend.margin = ggplot2::margin(t = -0.5, unit = "cm"),
        legend.title = ggplot2::element_blank()
      ) +
      ggplot2::guides(color = ggplot2::guide_legend(
        keywidth = 0.0,
        keyheight = 0.15,
        default.unit = "inch",
        override.aes = list(size = 1.25)
      ))
  }

  return(p_out)
}


make_volcano_plot <- function(discovery_result, p_thresh, x_limits = c(-1.5, 1.5), transparency = 0.5, point_size = 0.55) {
  p_lower_lim <- 1e-12
  out <- ggplot2::ggplot(
    data = discovery_result |> dplyr::mutate(
      reject = p_value < p_thresh,
      p_value = ifelse(p_value < p_lower_lim, p_lower_lim, p_value),
      log_2_fold_change = ifelse(log_2_fold_change > x_limits[2], x_limits[2], log_2_fold_change),
      log_2_fold_change = ifelse(log_2_fold_change < x_limits[1], x_limits[1], log_2_fold_change)
    ),
    mapping = ggplot2::aes(x = log_2_fold_change, y = p_value, col = reject)
  ) +
    ggplot2::geom_point(alpha = transparency, size = point_size) +
    ggplot2::scale_y_continuous(trans = revlog_trans(10), expand = c(0.02, 0)) +
    get_my_theme() +
    ggplot2::xlab("Log fold change") +
    ggplot2::ylab("P-value") +
    (if (!is.na(p_thresh)) ggplot2::geom_hline(yintercept = p_thresh, linetype = "dashed") else NULL) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_color_manual(values = c("dodgerblue", "blueviolet")) +
    ggplot2::ggtitle("Discovery volcano plot")
  return(out)
}


#' Plot run discovery analysis
#'
#' `plot_run_discovery_analysis()` creates a visualization of the outcome of the discovery analysis. The visualization consists of four plots:
#' -   The upper left plot superimposes the discovery p-values (blue) on top of the negative control p-values (red) on an untransformed scale.
#' -   The upper right plot is the same as the upper left plot, but the scale is negative log-10 transformed. The discovery p-values ideally should trend above the diagonal line, indicating the presence of signal in the discovery set. The horizontal dashed line indicates the multiple testing threshold; discovery pairs whose p-value falls above this line are called as significant.
#' -   The bottom left panel is a volcano plot of the p-values and log fold changes of the discovery pairs. Each point corresponds to a pair; the estimated log-2 fold change of the pair is plotted on the horizontal axis, and the (negative log-10 transformed) p-value is plotted on the vertical axis. The horizontal dashed line again indicates the multiple testing threshold. Points above the dashed line (colored in purple) are called as discoveries, while points below (colored in blue) are called as insignificant.
#' -   The bottom right panel is a text box displaying the number of discovery pairs called as significant.
#' @param sceptre_object a \code{sceptre_object} that has had \code{run_discovery_analysis} called on it
#' @param x_limits (optional; default \code{c(-1.5, 1.5)}) a numeric vector of length 2 giving the lower and upper limits of the x-axis (corresponding to log-2 fold change) for the "Discovery volcano plot" panel
#' @param point_size (optional; default \code{0.55}) the size of the individual points in the plot
#' @param transparency (optional; default \code{0.8}) the transparency of the individual points in the plot
#' @param return_indiv_plots (optional; default \code{FALSE}) if \code{FALSE} then a list of \code{ggplot} is returned; if \code{TRUE} then a single \code{cowplot} object is returned.
#' @return  a single \code{cowplot} object containing the combined panels (if \code{return_indiv_plots} is set to \code{TRUE}) or a list of the individual panels (if \code{return_indiv_plots} is set to \code{FALSE})
#'
#' @export
#' @examples
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' # import data
#' sceptre_object <- import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' )
#' positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
#' discovery_pairs <- construct_cis_pairs(sceptre_object,
#'   positive_control_pairs = positive_control_pairs,
#'   distance_threshold = 5e6
#' )
#' sceptre_object |>
#'   set_analysis_parameters(
#'     side = "left",
#'     discovery_pairs = discovery_pairs,
#'     resampling_mechanism = "permutations",
#'   ) |>
#'   assign_grnas(method = "thresholding") |>
#'   run_qc() |>
#'   run_calibration_check(
#'     parallel = TRUE,
#'     n_processors = 2
#'   ) |>
#'   run_discovery_analysis(
#'     parallel = TRUE,
#'     n_processors = 2
#'   ) |>
#'   plot_run_discovery_analysis()
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
  p1 <- compare_calibration_and_discovery_results(
    calibration_result = calibration_result,
    discovery_result = discovery_result,
    p_thresh = p_thresh,
    point_size = point_size,
    transparency = transparency,
    transform_scale = FALSE,
    include_legend = TRUE,
    include_y_axis_text = TRUE
  )
  # create tail qq plot
  p2 <- compare_calibration_and_discovery_results(
    calibration_result = calibration_result,
    discovery_result = discovery_result,
    p_thresh = p_thresh,
    point_size = point_size,
    transparency = transparency,
    transform_scale = TRUE,
    include_legend = FALSE,
    include_y_axis_text = FALSE
  )
  # make the volcano plot
  discovery_result_downsample <- downsample_result_data_frame(result_df = discovery_result)
  p3 <- make_volcano_plot(
    discovery_result = discovery_result_downsample,
    p_thresh = p_thresh,
    transparency = transparency,
    point_size = point_size,
    x_limits = x_limits
  )
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
#' `plot_covariates()` creates a histogram of the covariates `response_n_nonzero`, `response_n_umis`, and (if applicable) `response_p_mito`. Cellwise QC removes cells that lie in the extreme right tail of the `response_p_mito` distribution or that lie in the extreme left *or* right tail of the `response_n_nonzero` or `response_n_umis` distribution. To help guide the selection of QC thresholds, `plot_covariates()` plots candidate QC thresholds as vertical lines on the histograms. The optional arguments `response_n_nonzero_range`, `response_n_umis_range`, and `p_mito_threshold` control the location of these candidate QC thresholds. `response_n_nonzero_range` (resp., `response_n_umis_range`) is a length-two vector of quantiles (default: `c(0.01, 0.99)`) indicating the location at which to draw candidate QC thresholds on the `response_n_nonzero` (resp., `response_n_umis`) histogram. Next, `p_mito_threshold` is a single numeric value in the interval \[0,1\] specifying the location at which to draw a candidate QC threshold on the `response_p_mito` plot.
#'
#' @note If `run_qc()` has already been called on the `sceptre_object`, then the parameters `response_n_umis_range`, `response_n_nonzero_range`, and `p_mito_threshold` are set to the corresponding parameters within the `sceptre_object`.
#'
#' @param sceptre_object a `sceptre_object`
#' @param response_n_umis_range (optional; default \code{c(0.01, 0.99)}) a length-2 vector of quantiles indicating the location at which to draw vertical lines on the `response_n_umis` histogram
#' @param response_n_nonzero_range (optional; default \code{c(0.01, 0.99)}) a length-2 vector of quantiles indicating the location at which to draw vertical lines on the `response_n_nonzero` histogram
#' @param p_mito_threshold (optional; default \code{0.2}) a single numeric value in the interval \[0,1\] specifying the location at which to draw a vertical line on the `response_p_mito` histogram. Note that `p_mito_threshold` is an absolute number rather than a percentile.
#' @param return_indiv_plots (optional; default \code{FALSE}) if \code{FALSE}, then a list of \code{ggplot} objects is returned; if \code{TRUE} then a single \code{cowplot} object is returned.
#'
#' @return a single \code{cowplot} object containing the combined panels (if \code{return_indiv_plots} is set to \code{TRUE}) or a list of the individual panels (if \code{return_indiv_plots} is set to \code{FALSE})
#' @export
#' @examples
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' ) |> plot_covariates(p_mito_threshold = 0.07)
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
    p1 <- ggplot2::ggplot(
      data = data.frame(x = v),
      mapping = ggplot2::aes(x = x)
    ) +
      ggplot2::geom_histogram(col = "darkblue", fill = "grey90", bins = 50) +
      get_my_theme() +
      ggplot2::scale_y_continuous(expand = c(0, NA)) +
      ggplot2::geom_vline(xintercept = cutoffs[1], col = "darkorchid1", lwd = 1.0) +
      (if (length(cutoffs) == 2) {
        ggplot2::geom_vline(xintercept = cutoffs[2], col = "darkorchid1", lwd = 1.0)
      } else {
        NULL
      }) +
      ggplot2::ggtitle(plot_tit) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank())
  }
  p1 <- make_histogram(covariate_data_frame$response_n_nonzero, response_n_nonzero_range,
    "Response N nonzero",
    use_quantile = TRUE
  )
  p2 <- make_histogram(covariate_data_frame$response_n_umis, response_n_umis_range,
    "Response N UMIs",
    use_quantile = TRUE
  )
  p_mito_present <- "response_p_mito" %in% colnames(covariate_data_frame)
  if (p_mito_present) {
    p3 <- make_histogram(covariate_data_frame$response_p_mito, p_mito_threshold,
      plot_tit = "Percent mito", use_quantile = FALSE
    )
    p_out <- if (return_indiv_plots) list(p1, p2, p3) else cowplot::plot_grid(p1, p2, p3, NULL, ncol = 2)
  } else {
    p_out <- if (return_indiv_plots) list(p1, p2) else cowplot::plot_grid(p1, p2, ncol = 2)
  }
  return(p_out)
}


#' Plot run QC
#'
#' `plot_run_qc()` creates a visualization of the outcome of the QC step. The top panel depicts the outcome of the cellwise QC. The various cellwise QC filters (e.g., "N nonzero responses," "N response UMIs," "Percent mito", etc.) are shown on the horizontal axis, and the percentage of cells removed due application of a given QC filter is shown on the vertical axis. Note that a cell can be flagged by multiple QC filters; for example, a cell might have an extremely high `response_n_umi` value *and* an extremely high `response_n_nonzero` value. Thus, the height of the "any filter" bar (which indicates the percentage of cells removed due to application of *any* filter) need not be equal to the sum of the heights of the other bars. The bottom panel depicts the outcome of the pairwise QC. Each point corresponds to a target-response pair; the vertical axis (resp., horizontal axis) indicates the `n_nonzero_trt` (resp., `n_nonzero_cntrl`) value of that pair. Pairs for which `n_nonzero_trt` or `n_nonzero_cntrl` fall below the threshold are removed (red), while the remaining pairs are retained (green).
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
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' # import data
#' sceptre_object <- import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' )
#' discovery_pairs <- construct_cis_pairs(sceptre_object)
#' sceptre_object |>
#'   set_analysis_parameters(
#'     discovery_pairs = discovery_pairs,
#'     side = "left"
#'   ) |>
#'   assign_grnas(method = "thresholding") |>
#'   run_qc() |>
#'   plot_run_qc()
plot_run_qc <- function(sceptre_object, downsample_pairs = 10000L, point_size = 0.55, transparency = 0.8, return_indiv_plots = FALSE) {
  if (!sceptre_object@functs_called["run_qc"]) {
    stop("This `sceptre_object` has not yet had `run_qc` called on it.")
  }
  p_a <- plot_cellwise_qc(sceptre_object)
  if (nrow(sceptre_object@discovery_pairs_with_info) >= 1L) {
    p_b <- plot_pairwise_qc(sceptre_object, downsample_pairs = 10000L, point_size = 0.55, transparency = 0.8, return_indiv_plots = FALSE)
    # combine the plots
    if (return_indiv_plots) {
      p_out <- list(p_a, p_b)
    } else {
      p_out <- cowplot::plot_grid(p_a, p_b, ncol = 1)
    }
  } else {
    p_out <- p_a
  }
  return(p_out)
}

plot_cellwise_qc <- function(sceptre_object) {
  cell_removal_metrics <- sceptre_object@cell_removal_metrics
  n_orig_cells <- ncol(get_response_matrix(sceptre_object))
  cell_removal_metrics_frac <- cell_removal_metrics / n_orig_cells * 100
  df <- data.frame(
    fraction_cells_removed = cell_removal_metrics_frac,
    Filter = c(
      "N response UMIs", "N nonzero responses", "Percent mito",
      "Zero or 2+ gRNAs", "User-specified", "Any filter"
    )
  )
  if (!sceptre_object@low_moi) df <- df |> dplyr::filter(Filter != "Zero or 2+ gRNAs")
  # make a barplot. remove x-axis text
  p_a <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = Filter, y = fraction_cells_removed)) +
    ggplot2::geom_bar(stat = "identity", fill = "grey90", col = "darkblue") +
    get_my_theme() +
    ggplot2::scale_y_continuous(expand = c(0, NA)) +
    ggplot2::ylab("Percent cells removed") +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 20)) +
    get_my_theme() +
    ggplot2::theme(legend.position = "none", axis.title.x = ggplot2::element_blank()) +
    ggplot2::ggtitle("Cellwise QC")
  return(p_a)
}

plot_pairwise_qc <- function(sceptre_object, downsample_pairs = 10000L, point_size = 0.55, transparency = 0.8, return_indiv_plots = FALSE) {
  my_cols <- c("mediumseagreen", "indianred2")
  my_breaks <- c(0, 1, 3, 7, 50, 500, 5000, 50000)
  discovery_pairs <- sceptre_object@discovery_pairs_with_info |>
    dplyr::mutate(pass_qc = ifelse(pass_qc, "Pass", "Fail")) |>
    dplyr::mutate(pass_qc = factor(pass_qc, levels = c("Pass", "Fail")))
  if (nrow(discovery_pairs) >= downsample_pairs) discovery_pairs <- discovery_pairs |> dplyr::sample_n(downsample_pairs)

  p_b <- ggplot2::ggplot(
    data = discovery_pairs,
    mapping = ggplot2::aes(x = n_nonzero_trt, y = n_nonzero_cntrl, col = pass_qc)
  ) +
    ggplot2::geom_point(alpha = transparency, size = point_size) +
    get_my_theme() +
    ggplot2::scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10, sigma = 1),
      breaks = my_breaks
    ) +
    ggplot2::scale_x_continuous(
      trans = scales::pseudo_log_trans(base = 10, sigma = 1),
      breaks = my_breaks
    ) +
    ggplot2::geom_hline(yintercept = sceptre_object@n_nonzero_cntrl_thresh) +
    ggplot2::geom_vline(xintercept = sceptre_object@n_nonzero_trt_thresh) +
    ggplot2::xlab("N nonzero trt. cells") +
    ggplot2::ylab("N nonzero cntrl. cells") +
    ggplot2::theme(legend.position = "bot") +
    ggplot2::scale_color_manual(values = my_cols) +
    ggplot2::ggtitle("Pairwise QC") +
    ggplot2::theme(
      legend.position.inside = c(0.85, 0.2),
      legend.margin = ggplot2::margin(t = -0.5, unit = "cm"),
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(
      keywidth = 0.0,
      keyheight = 0.15,
      default.unit = "inch",
      override.aes = list(size = 1.25)
    ))
  return(p_b)
}

###############
# 6. PLOT POWER
###############
#' Plot run power check
#'
#' \code{plot_run_power_check()} creates a visualization of the outcome of the power check analysis. Each point in the plot corresponds to a target-response pair, with positive control pairs in the left column and negative control pairs in the right column. The vertical axis indicates the p-value of a given pair; smaller (i.e., more significant) p-values are positioned higher along this axis (p-values truncated at `clip_to` for visualization). The positive control p-values should be small, and in particular, smaller than the negative control p-values.
#'
#' @param sceptre_object a \code{sceptre_object} that has had \code{run_power_check()} called on it
#' @param point_size (optional; default \code{1}) the size of the individual points in the plot
#' @param transparency (optional; default \code{0.8}) the transparency of the individual points in the plot
#' @param clip_to (optional; default \code{1e-20}) p-values smaller than this value are set to \code{clip_to} for better visualization. If \code{clip_to=0} is used then no clipping is done.
#'
#' @return a single \code{ggplot2} plot
#' @export
#' @examples
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' # import data
#' sceptre_object <- import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' )
#' positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
#' sceptre_object |>
#'   set_analysis_parameters(
#'     side = "left",
#'     positive_control_pairs = positive_control_pairs,
#'     resampling_mechanism = "permutations",
#'   ) |>
#'   assign_grnas(method = "thresholding") |>
#'   run_qc() |>
#'   run_calibration_check(
#'     parallel = TRUE,
#'     n_processors = 2,
#'     n_calibration_pairs = 500,
#'     calibration_group_size = 2
#'   ) |>
#'   run_power_check() |>
#'   plot_run_power_check()
plot_run_power_check <- function(sceptre_object, point_size = 1, transparency = 0.8, clip_to = 1e-20) {
  if (!sceptre_object@functs_called["run_power_check"]) {
    stop("This `sceptre_object` has not yet had `run_power_check` called on it.")
  }
  my_theme <- get_my_theme()
  set.seed(3)
  my_cols <- c("mediumseagreen", "firebrick1")

  pos_ctrl_pvals <- sceptre_object@power_result$p_value |> stats::na.omit()
  neg_ctrl_pval_sub <- if (nrow(sceptre_object@calibration_result) >= 1) {
    downsample_result_data_frame(
      result_df = sceptre_object@calibration_result
    ) |> dplyr::pull(p_value)
  } else {
    numeric(0)
  }
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
      title = "Positive and negative control p-values"
    ) +
    my_theme

  return(p)
}


downsample_result_data_frame <- function(result_df, downsample_pairs = 1000) {
  result_df |>
    dplyr::mutate(bin = cut(p_value, breaks = c(10^(-seq(0, 20)), 0))) |>
    dplyr::group_by(bin) |>
    dplyr::sample_n(size = min(dplyr::n(), downsample_pairs))
}


#' Plot response-gRNA-target pair
#'
#' `plot_response_grna_target_pair()` creates a violin plot of the expression level of a given response as a function of the "treatment status" (i.e., treatment or control) of a given gRNA target. The left (resp., right) violin plot shows the expression level of the response in treatment (resp., control) cells. The expression level is normalized by dividing by `n_response_umis`, adding a pseudo-count of 1 and then taking the log transform. If the given response-gRNA-target pair has been analyzed, the p-value for the test of association also is displayed.
#'
#' If `grna_integration_strategy` is set to `"singleton"`, then `grna_target` should be set to a gRNA ID.
#' @param sceptre_object a `sceptre_object` that has had `run_qc()` called on it
#' @param response_id a string containing a response ID
#' @param grna_target a string containing a gRNA target (or, if `grna_integration_strategy` is set to `"singleton"`, an individual gRNA ID)
#'
#' @return a violin plot
#' @export
#' @examples
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' # import data
#' sceptre_object <- import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' )
#' discovery_pairs <- construct_cis_pairs(sceptre_object)
#' sceptre_object |>
#'   set_analysis_parameters(
#'     side = "left",
#'     discovery_pairs = discovery_pairs,
#'     resampling_mechanism = "permutations",
#'   ) |>
#'   assign_grnas(method = "thresholding") |>
#'   run_qc() |>
#'   run_discovery_analysis(
#'     parallel = TRUE,
#'     n_processors = 2
#'   ) |>
#'   plot_response_grna_target_pair(
#'     response_id = "ENSG00000211641",
#'     grna_target = "candidate_enh_15"
#'   )
plot_response_grna_target_pair <- function(sceptre_object, response_id, grna_target) {
  # check that grnas have been assigned and qc has been called
  functs_called <- sceptre_object@functs_called
  singleton_integration_strategy <- sceptre_object@grna_integration_strategy == "singleton"
  if (!functs_called[["assign_grnas"]]) {
    stop("`assign_grnas()` must be called on the `sceptre_object`.")
  }
  if (!functs_called[["run_qc"]]) {
    stop("`run_qc()` must be called on the `sceptre_object`.")
  }

  # get grna integration strategy, control group
  control_group_complement <- sceptre_object@control_group_complement
  cells_in_use <- sceptre_object@cells_in_use

  # check that grna_target is present in grna_assignments
  grna_group_idxs <- sceptre_object@grna_assignments$grna_group_idxs
  grna_group_names <- names(grna_group_idxs)
  if (!grna_target %in% grna_group_names) {
    stop(if (singleton_integration_strategy) "gRNA ID" else "gRNA target `", grna_target, "` is not present within the data.")
  }

  # check that the response is present within the data
  response_matrix <- get_response_matrix(sceptre_object)
  if (!(response_id %in% rownames(response_matrix))) {
    stop("Response `", response_id, "` is not present within the data.")
  }

  # extract the counts for this pair; filter for cells in use
  y <- load_row(mat = response_matrix, id = response_id)[cells_in_use]

  # get the sequencing depths; filter for cells in use
  response_n_umis <- sceptre_object@covariate_data_frame$response_n_umis[cells_in_use]

  # compute the normalized counts
  normalized_counts <- log(10000 * y / response_n_umis + 1)

  # get the treated cells and control cells
  trt_cells <- normalized_counts[grna_group_idxs[[grna_target]]]
  if (control_group_complement) { # complement set
    cntrl_cells <- normalized_counts[-grna_group_idxs[[grna_target]]]
  } else { # nt cells
    nt_idxs <- sceptre_object@grna_assignments$all_nt_idxs
    cntrl_cells <- normalized_counts[nt_idxs]
  }

  # construct the data frame to plot
  to_plot <- data.frame(
    normalized_count = c(trt_cells, cntrl_cells),
    treatment = c(
      rep("Treatment", length(trt_cells)),
      rep("Control", length(cntrl_cells))
    ) |>
      factor(levels = c("Treatment", "Control"))
  )

  # obtain the p-value (if available); deal with the singleton situation
  p_val <- 1.5
  df_list <- list(if (functs_called[["run_power_check"]]) sceptre_object@power_result else NULL,
                  if (functs_called[["run_discovery_analysis"]]) sceptre_object@discovery_result else NULL) |>
    purrr::compact()
  for (curr_df in df_list) {
    if (singleton_integration_strategy) {
      subset_result <- curr_df |>
        dplyr::filter(response_id == !!response_id, grna_id == !!grna_target)
    } else {
      subset_result <- curr_df |>
        dplyr::filter(response_id == !!response_id, grna_target == !!grna_target)
    }
    if (nrow(subset_result) >= 1) {
      p_val <- subset_result$p_value
      if (is.na(p_val)) p_val <- 1
    }
  }

  # obtain the annotation
  annotation <- cut(p_val,
    breaks = c(2, 1, 10^(-seq(2, 10)), 0),
    labels = c(paste0("p <= 10^{", seq(-10, -2), "}"), "p > 0.01", "")
  ) |> as.character()
  # create the plot
  set.seed(4)
  to_plot_downsample <- to_plot |>
    dplyr::mutate(is_zero = (normalized_count == 0)) |>
    dplyr::group_by(is_zero, treatment) |>
    dplyr::sample_n(size = min(dplyr::n(), 1000))
  p_out <- ggplot2::ggplot(data = to_plot, mapping = ggplot2::aes(x = treatment, y = normalized_count, col = treatment)) +
    ggplot2::geom_violin(draw_quantiles = 0.5, linewidth = 0.6) +
    ggplot2::geom_jitter(data = to_plot_downsample, alpha = 0.1, size = 0.5) +
    get_my_theme() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_color_manual(values = c("dodgerblue4", "firebrick4")) +
    ggplot2::xlab("Treatment status") +
    ggplot2::ylab("Normalized expression") +
    ggplot2::annotate("text", x = 1.5, y = max(to_plot$normalized_count) + 0.5, label = annotation, parse = TRUE) +
    ggplot2::scale_y_continuous(expand = c(0.0, 0.1), limits = c(-0.01, max(to_plot$normalized_count) + 0.7)) +
    ggplot2::ggtitle(paste0("Response: ", response_id, "\ngRNA", if (singleton_integration_strategy) "" else " target", ": ", grna_target)) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10))

  return(p_out)
}
