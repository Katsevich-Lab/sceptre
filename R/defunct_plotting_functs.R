make_n_nonzero_cntrl_vs_trt_cells_plot <- function(sceptre_object) {
  my_cols <- c("mediumseagreen", "indianred2")
  my_breaks <- c(0, 1, 3, 7, 50, 500, 5000, 50000)
  discovery_pairs <- sceptre_object@discovery_pairs_with_info |>
    dplyr::mutate(pass_qc = ifelse(pass_qc, "Pass", "Fail")) |>
    dplyr::mutate(pass_qc = factor(pass_qc, levels = c("Pass", "Fail")))
  ggplot2::ggplot(data = discovery_pairs,
                  mapping = ggplot2::aes(x = n_nonzero_trt, y = n_nonzero_cntrl, col = pass_qc)) +
    ggplot2::geom_point(alpha = 0.8, size = 0.8) + get_my_theme() +
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
    ggplot2::ggtitle("Pairwise QC plot") +
    ggplot2::theme(legend.position = c(0.85, 0.15),
                   legend.margin = ggplot2::margin(t = -0.5, unit = "cm"),
                   legend.title = ggplot2::element_blank()) +
    ggplot2::guides(color = ggplot2::guide_legend(
      keywidth = 0.0,
      keyheight = 0.15,
      default.unit = "inch",
      override.aes = list(size = 1.25)))
}

make_cells_per_grna_group_plot <- function(sceptre_object) {
  n_cells_per_grna_group <- sapply(sceptre_object@grna_assignments$grna_group_idxs, length)
  df <- data.frame(x = n_cells_per_grna_group,
                   y = names(sceptre_object@grna_assignments$grna_group_idxs)) |>
    dplyr::arrange(n_cells_per_grna_group) |>
    dplyr::mutate(y = factor(y, labels = y, levels = y))
  p <- ggplot2::ggplot(data = df,
                       mapping = ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_bar(stat = "identity", width = 0.5, fill = "darkblue", col = "black") +
    get_my_theme() + ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::xlab("N cells") + ggplot2::ylab("gRNA group") +
    ggplot2::ggtitle("N cells per gRNA group")
  return(list(plot = p, median = median(n_cells_per_grna_group)))
}

make_n_grna_groups_per_cell_plot <- function(sceptre_object) {
  n_grna_groups_per_cell <- sceptre_object@grna_assignments$grna_group_idxs |>
    unlist() |> tabulate()
  df <- data.frame(x = n_grna_groups_per_cell)
  p <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = x)) +
    ggplot2::geom_histogram(binwidth = 1, fill = "grey90", color = "darkblue") +
    ggplot2::scale_y_continuous(trans = "log10", expand = c(0, 0)) +
    get_my_theme() + ggplot2::ylab("Count") +
    ggplot2::ggtitle("N gRNA groups per cell") +
    ggplot2::xlab("gRNA groups per cell")
  return(list(plot = p, median = median(n_grna_groups_per_cell)))
}

plot_prepare_analysis <- function(sceptre_object, return_indiv_plots = FALSE) {
  # perform checking
  if (!sceptre_object@analysis_prepared) {
    stop("Analysis must be prepared before calling `plot_prepare_analysis` function.")
  }

  # make plot
  p_a_out <- make_cells_per_grna_group_plot(sceptre_object)
  str <- paste0("Median cells per\ngRNA group: ", p_a_out$median)
  p_a <- p_a_out$plot
  if (!sceptre_object@low_moi) {
    p_b_out <- make_n_grna_groups_per_cell_plot(sceptre_object)
    str <- paste0(str, "\n\nMedian gRNA groups\nper cell: ", p_b_out$median)
    p_b <- p_b_out$plot
  }
  p_c <- ggplot2::ggplot() +
    ggplot2::annotate(geom = "text", label = str, x = 1.0, y = 1.0) +
    ggplot2::theme_void() +
    ggplot2::xlim(c(0, 2)) +
    ggplot2::ylim(c(0, 2))
  p_d <- make_n_nonzero_cntrl_vs_trt_cells_plot(sceptre_object)

  if (return_indiv_plots) {
    p_out <- if (sceptre_object@low_moi) list(p_a, p_c, p_d) else list(p_a, p_b, p_c, p_d)
  } else {
    if (sceptre_object@low_moi) {
      p_out <- cowplot::plot_grid(cowplot::plot_grid(p_a, p_c, nrow = 1, rel_widths = c(0.65, 0.35)), p_d, nrow = 2)
    } else {
      p_out <- cowplot::plot_grid(p_a, p_b, p_c, p_d, nrow = 2, align = "h")
    }
  }
  return(p_out)
}


plot_covariates <- function(sceptre_object, return_indiv_plots = FALSE) {
  # obtain the covariate df, user covariates, and auto covariates
  covariate_df <- sceptre_object@covariate_data_frame
  user_specified_covariates <- sceptre_object@user_specified_covariates

  # get the user-specified and automatic covariates to plot
  set.seed(4L)
  user_covariates_to_plot <- sample(x = user_specified_covariates, replace = FALSE,
                                    size = min(2L, length(user_specified_covariates)))
  auto_covariate_ordering <- c("response_n_umis", "grna_n_umis",
                               if ("response_p_mito" %in% colnames(covariate_df)) "response_p_mito" else NULL,
                               "response_n_nonzero", "grna_n_nonzero")
  auto_covariates_to_plot <- auto_covariate_ordering[seq(1L, 4L - length(user_covariates_to_plot))]

  # histogram plotting code
  make_histogram <- function(vect, log_trans_x, plot_name) {
    zero_present <- any(vect == 0)
    p <- data.frame(vect = vect) |>
      ggplot2::ggplot(ggplot2::aes(x = vect)) +
      ggplot2::geom_histogram(bins = 30, color = "darkblue", fill = "grey90") +
      get_my_theme() +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::ggtitle(plot_name)
    if (log_trans_x) {
      if (zero_present) {
        p <- p + ggplot2::scale_x_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 10),
                                             n.breaks = 4L)
      } else {
        p <- p +  ggplot2::scale_x_continuous(trans = "log10")
      }
    }

    return(p)
  }

  # barplot plotting code
  make_barplot <- function(vect, plot_name) {
    p <- data.frame(vect = vect) |>
      ggplot2::ggplot(ggplot2::aes(x = vect, fill = vect)) +
      ggplot2::geom_bar(col = "black") +
      get_my_theme() +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank()) +
      ggplot2::ggtitle(plot_name)
    return(p)
  }

  # combine the covariates
  covariates_to_plot <- c(user_covariates_to_plot, auto_covariates_to_plot)
  ps <- sapply(covariates_to_plot, FUN = function(curr_covariate) {
    vect <- covariate_df[,curr_covariate]
    tit <- gsub(x = curr_covariate, pattern = "_", fixed = TRUE, replacement = " ")
    if (is(vect, "numeric")) {
      make_histogram(vect, grepl(pattern = "n_umis|n_nonzero", x = curr_covariate), tit)
    } else {
      make_barplot(vect, tit)
    }
  }, simplify = FALSE)

  if (return_indiv_plots) {
    out <- ps
  } else {
    out <- cowplot::plot_grid(plotlist = ps, nrow = 2, labels = "auto", align = "h")
  }
  return(out)
}
