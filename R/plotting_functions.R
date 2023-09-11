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
#' @export
plot_grna_count_distributions <- function(sceptre_object, n_grnas_to_plot = 4L, grnas_to_plot = NULL) {
  grna_matrix <- sceptre_object@grna_matrix
  if (is.null(grnas_to_plot)) {
    grnas_to_plot <- sample(x = rownames(grna_matrix), size = min(nrow(grna_matrix), n_grnas_to_plot), replace = FALSE)
  }
  grna_matrix <- set_matrix_accessibility(grna_matrix, make_row_accessible = TRUE)
  grna_expressions <- lapply(X = grnas_to_plot, function(grna_id) {
    load_csr_row(j = grna_matrix@j, p = grna_matrix@p, x = grna_matrix@x,
                 row_idx = which(grna_id == rownames(grna_matrix)), n_cells = ncol(grna_matrix))
  }) |> unlist()
  grna_ids_rep <- rep(factor(grnas_to_plot), each = ncol(grna_matrix))
  to_plot <- data.frame(grna_id = grna_ids_rep, grna_expressions = grna_expressions) |>
    dplyr::filter(grna_expressions < 10000)
  p <- ggplot2::ggplot(data = to_plot, mapping = ggplot2::aes(x = grna_expressions)) +
    ggplot2::geom_histogram(bins = 20, col = "midnightblue", fill = "grey90") + get_my_theme() +
    ggplot2::facet_wrap(grna_id ~ ., scales = "free_x", nrow = floor(sqrt(length(grnas_to_plot)))) +
    ggplot2::scale_x_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 1),
                                breaks = c(0, 1, 3, 7, 50, 500)) +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 0.5),
                                breaks = c(0, 10, 100, 1000, 100000), expand = c(0, NA)) +
    ggplot2::xlab("gRNA count") + ggplot2::ylab("N cells")
  return(p)
}

#' Plot the number of gRNAs per cell
#'
#' @param n_grnas_to_plot (optional; default \code{2}) the number of different gRNAs shown in the plots of gRNA count versus cell assigment.
#' @param grnas_to_plot (optional; default \code{NULL}) the names of specific gRNAs to plot; if \code{NULL} then gRNAs are chosen at random.
#' @param return_indiv_plots (optional; default \code{FALSE}) return a single combined plot (\code{TRUE}) or the individual plots in a list (\code{FALSE})?
#' @param n_max_0_grna_unprtb_plot (optional; default \code{1000}) there may be many cells with a gRNA count of 0 in the unperturbed group. This can slow down the plotting without adding useful information, so at most \code{n_max_0_grna_unprtb_plot} points from this group are plotted. Setting this to \code{Inf} will guarantee no downsampling occurs.
#'
#' @return a single \code{cowplot} object containing the combined panels (if \code{return_indiv_plots} is set to \code{TRUE}) or a list of the individual panels (if \code{return_indiv_plots} is set to \code{FALSE}).
#' @export
plot_assign_grnas <- function(sceptre_object, n_grnas_to_plot = 2L, grnas_to_plot = NULL, transparency = 0.8, point_size = 0.9, return_indiv_plots = FALSE, n_max_0_grna_unprtb_plot = 1000) {
  init_assignments <- sceptre_object@initial_grna_assignment_list
  grna_matrix <- sceptre_object@grna_matrix |> set_matrix_accessibility(make_row_accessible = TRUE)
  grna_ids <- rownames(grna_matrix)
  lowmoi <- sceptre_object@low_moi
  # sample grnas to plot
  if (is.null(grnas_to_plot)) grnas_to_plot <- sample(grna_ids, size = min(n_grnas_to_plot, length(grna_ids)), replace = FALSE)

  to_plot_a <- lapply(X = grnas_to_plot, function(grna_id) {
    assignment <- multiple_grnas <- logical(length = ncol(grna_matrix)) # logical vecs w/ one entry per cell
    assignment[init_assignments[[grna_id]]] <- TRUE # for this grna, `assignment` indicates which cells got this grna initially
    multiple_grnas[sceptre_object@cells_w_multiple_grnas] <- TRUE  # indicates which cells have >1 grna
    g <- grna_matrix[grna_id,] # all cell UMI counts for this grna
    data.frame(g = g,
               assignment = ifelse(assignment, "perturbed", "unperturbed") |> factor(),
               grna_id = grna_id |> factor(),
               # making this a factor so the ordering never changes
               multiple_grnas = ifelse(multiple_grnas, "yes", "no") |> factor(levels=c("no","yes")))
  }) |> data.table::rbindlist()

  # downsampling the unperturbed cells with 0 grna expression, if the number of those cells
  # exceeds `n_max_0_grna_unprtb_plot`.
  is_0_grna_and_unperturbed <- with(to_plot_a, assignment == "unperturbed" & g == 0)
  if(sum(is_0_grna_and_unperturbed) > n_max_0_grna_unprtb_plot) {
    idx_to_remove <- rownames(to_plot_a)[is_0_grna_and_unperturbed] |>
      sample(sum(is_0_grna_and_unperturbed) - n_max_0_grna_unprtb_plot)
    to_plot_a <- to_plot_a[! rownames(to_plot_a) %in% idx_to_remove, ]
  }

  # plot a
  p_a <- ggplot2::ggplot(data = to_plot_a, mapping = (if (lowmoi) {
    ggplot2::aes(x = assignment, y = g, col = assignment, shape = multiple_grnas)
  } else ggplot2::aes(x = assignment, y = g, col = assignment))) +
    ggplot2::geom_jitter(alpha = transparency, size = point_size) +
    ggplot2::facet_wrap(nrow = 1, facets = grna_id ~ ., scales = "free_y") +
    ggplot2::scale_y_continuous(trans = "log1p", breaks = c(0, 1, 5, 15, 50, 200, 1000, 5000, 10000)) +
    ggplot2::xlab("Assignment") +
    ggplot2::ylab("gRNA count") +
    get_my_theme() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
    ggplot2::scale_color_manual(values = c("firebrick1", "darkorchid1"))

  # in low moi we want to have a legend for the shape parameter
  # to save space, we'll do it inside the first panel and just for
  # one level of the variable
  if(lowmoi) {
    # setting up the data for making the legend as a geom
    annotate_df <- data.frame(
      assignment = "perturbed",
      g = .5,
      # making sure it's in the first panel
      grna_id = factor(levels(to_plot_a$grna_id)[1], levels = levels(to_plot_a$grna_id)),
      multiple_grnas = NA # just need this column here, but the color is overwritten
    )
    p_a <- p_a + ggplot2::geom_label(
      data = annotate_df,
      label = c("\U25B2 cells with\nmultiple gRNAs"),
      color="black",
      size = 2.5,
      label.padding = ggplot2::unit(0.1, "lines"),
    )
  }

  # plot b
  n_cells_per_grna <- sceptre_object@cells_per_grna
  to_plot_b <- data.frame(x = n_cells_per_grna,
                          y = names(n_cells_per_grna)) |>
    dplyr::arrange(n_cells_per_grna) |>
    dplyr::mutate(y = factor(y, labels = y, levels = y))
  p_b <- ggplot2::ggplot(data = to_plot_b,
                         mapping = ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_bar(stat = "identity", width = 0.5, fill = "grey90", col = "darkblue") +
    get_my_theme() + ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::xlab("N cells") + ggplot2::ylab("gRNA") +
    ggplot2::ggtitle("N cells per gRNA")

  # plot c
  if (sceptre_object@grna_assignment_method == "maximum") {
    p_c <- ggplot2::ggplot() +
      ggplot2::theme_void()
  } else {
    n_grnas_per_cell <- sceptre_object@grnas_per_cell
    moi <- mean(n_grnas_per_cell)
    to_plot_c <- data.frame(x = n_grnas_per_cell)
    p_c <- ggplot2::ggplot(data = to_plot_c, mapping = ggplot2::aes(x = x)) +
      ggplot2::geom_histogram(binwidth = 1, fill = "grey90", color = "darkblue") +
      ggplot2::scale_y_continuous(trans = "log1p",
                                  expand = c(0, 0),
                                  breaks = 10^(1:8)) +
      get_my_theme() + ggplot2::ylab("Frequency") +
      ggplot2::ggtitle("gRNAs per cell") +
      ggplot2::xlab("gRNAs per cell") +
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
  out
  return(out)
}

###########################
# 3. PLOT CALIBRATION CHECK
###########################
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
plot_run_calibration_check <- function(sceptre_object, return_indiv_plots = FALSE, point_size = 0.55, transparency = 0.8) {
  calibration_result <- sceptre_object@calibration_result
  if (nrow(calibration_result) == 0L) stop("Calibration check not yet called.")
  my_theme <- get_my_theme()

  # compute n rejections
  n_rejections <- sum(calibration_result$significant)
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
    ggplot2::geom_histogram(binwidth = 0.02, fill = "grey90", col = "black", boundary = 0) +
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
    ggplot2::scale_y_continuous(trans = revlog_trans(10), limits = c(1, 1e-12), expand = c(0.02, 0)) +
    get_my_theme() + ggplot2::xlab("Log fold change") + ggplot2::ylab("P-value") +
    (if (!is.na(p_thresh)) ggplot2::geom_hline(yintercept = p_thresh, linetype = "dashed") else NULL) +
    ggplot2::theme(legend.position = "none") + ggplot2::scale_color_manual(values = c("dodgerblue", "blueviolet")) +
    ggplot2::ggtitle("Discovery volcano plot")
  return(out)
}


#' @export
plot_run_discovery_analysis <- function(sceptre_object, return_indiv_plots = FALSE, x_limits = c(-1.5, 1.5), transparency = 0.8, point_size = 0.55) {
  # first, compute the rejection set
  if (nrow(sceptre_object@discovery_result) == 0L) stop("Discovery analysis not run.")
  discovery_result <- sceptre_object@discovery_result |> na.omit()
  calibration_result <- sceptre_object@calibration_result
  discovery_set <- discovery_result |> dplyr::filter(significant)
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
  str <- paste0("Number of discovery pairs\ncalled as significant\n(at alpha ", signif(sceptre_object@multiple_testing_alpha, 1), "): ", n_rejections, " of ", n_pairs)
  p4 <- ggplot2::ggplot() +
    ggplot2::annotate(geom = "text", label = str, x = 1.1, y = 1.2) +
    ggplot2::theme_void() +
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
#' @export
plot_run_qc <- function(sceptre_object, return_indiv_plots = FALSE, transparency = 0.8, point_size = 0.55) {
  my_cols <- c("mediumseagreen", "indianred2")
  my_breaks <- c(0, 1, 3, 7, 50, 500, 5000, 50000)
  discovery_pairs <- sceptre_object@discovery_pairs_with_info |>
    dplyr::mutate(pass_qc = ifelse(pass_qc, "Pass", "Fail")) |>
    dplyr::mutate(pass_qc = factor(pass_qc, levels = c("Pass", "Fail")))

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
  n_orig_cells <- ncol(sceptre_object@response_matrix)
  cell_removal_metrics_frac <- cell_removal_metrics/n_orig_cells * 100
  df <- data.frame(fraction_cells_removed = cell_removal_metrics_frac,
                   Filter = c("N response UMIs", "N nonzero responses", "percent mito", "multiple gRNAs", "user-specified", "any filter"))
  if (!sceptre_object@low_moi) df <- df |> dplyr::filter(Filter != "multiple gRNAs")
  # make a barplot. remove x-axis text
  p_a <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = Filter, y = fraction_cells_removed)) +
    ggplot2::geom_bar(stat = "identity", fill = "grey90", col = "darkblue") + get_my_theme() +
    ggplot2::scale_y_continuous(expand = c(0, NA)) +
    ggplot2::ylab("Percent. cells removed") +
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
#' Plot power result
#'
#' The \code{plot_run_power_check} function helps to visualize the results of a power check.
#'
#' The plot shows positive and negative control p-values, plotted with jitter so each p-value is visible, on a reversed log10 scale.
#'
#' @param return_indiv_plots (optional; default \code{FALSE}) ignored; kept for compatibility with other plotting functions.
#' @param clip_to (optional; default \code{1e-20}) p-values smaller than this value are set to \code{clip_to}, for better visualization. If \code{clip_to=0} is used then no clipping is done.
#' @return a single \code{ggplot2} plot.
#' @export
plot_run_power_check <- function(sceptre_object, return_indiv_plots = FALSE, point_size = 1,
                                 transparency = 0.8,  clip_to = 1e-20) {
  calibration_result <- sceptre_object@power_result
  if (nrow(calibration_result) == 0L) stop("Power check not yet called.")
  my_theme <- get_my_theme()
  set.seed(3)
  my_cols <- c("mediumseagreen", "firebrick1")

  pos_ctrl_pvals <- sceptre_object@power_result$p_value
  neg_ctrl_pvals <- sceptre_object@calibration_result$p_value

  group_names <- c("Positive Control", "Negative Control")
  df <- data.frame(
    lab = rep(
      group_names,
      c(length(pos_ctrl_pvals), length(neg_ctrl_pvals))
    ) |>
      factor(levels = group_names),
    p_values = c(pos_ctrl_pvals, neg_ctrl_pvals) |>
      pmax(clip_to)
  )

  p <- ggplot2::ggplot(
    data = df,
    mapping = ggplot2::aes(x = lab, y = p_values, color=lab)
  )
  if(clip_to > 0) {
    p <- p + ggplot2::geom_hline(yintercept = clip_to, linetype = "dashed", color="grey90")
  }
  p <- p +
    ggplot2::geom_jitter(width = .25, height = 0, size = point_size, alpha = transparency) +
    ggplot2::scale_y_continuous(trans = revlog_trans(base = 10), expand = c(0.01, 0)) +
    ggplot2::scale_color_manual(values = my_cols, guide = "none") +
    ggplot2::labs(
      x = "Pair type",
      y = "p-value",
      title = "Positive and negative control p-values") +
    my_theme

  return(p)
}
