# plotting functions

# taken from Julius Vainora's answer at https://stackoverflow.com/questions/54051576
# This is used to get information from a cowplot::plot_grid object
get_elements_from_plot_grid <- function(p, what) {
  unlist(sapply(p$layers, function(x) {
    idx <- which(x$geom_params$grob$layout$name == what)
    x$geom_params$grob$grobs[[idx]]$children[[1]]$label
  }))
}

# this test makes a single sceptre object and runs thru each of the main plots
# making sure their basic components are all there
# This is all done as one test just to save a little time on repeatedly creating the
# sceptre object, and since these aren't very intense tests
test_that("test all plots", {
  grna_target_data_frame <- data.frame(
    grna_id = c("id1", "id2", "id3", "nt1", "nt2"),
    grna_target = c("t1", "t2", "t3", "non-targeting", "non-targeting"),
    chr = "", start = 0, end = 1
  )
  num_grna <- nrow(grna_target_data_frame)
  num_cells <- 50
  num_responses <- 10

  set.seed(1)
  # using sample(0:1) so no entries can accidentally cross the threshold
  grna_matrix <- matrix(sample(0:1, num_grna * num_cells, replace = TRUE), num_grna, num_cells) |>
    `rownames<-`(grna_target_data_frame$grna_id)
  cells_expressing_t1 <- 1:10
  cells_expressing_t2 <- 11:20
  cells_expressing_t3 <- 21:30
  cells_expressing_nt1 <- 31:40
  cells_expressing_nt2 <- 41:50
  all_cells <- 1:num_cells

  grna_matrix["id1", cells_expressing_t1] <- 50
  grna_matrix["id2", cells_expressing_t2] <- 50
  # grna_matrix["id3", cells_expressing_t3] <- 50
  grna_matrix["nt1", cells_expressing_nt1] <- 50
  grna_matrix["nt2", cells_expressing_nt2] <- 50

  response_matrix <- matrix(rpois(num_responses * num_cells, 1), num_responses, num_cells) |>
    `rownames<-`(c("t1", "t2", "t3", paste0("response_", 4:num_responses)))

  response_matrix["t1", cells_expressing_t1] <- 100 # should be highly significant

  response_matrix["t2", ] <- 100 # should not be significant at all

  positive_control_pairs <- data.frame(
    grna_target = c("t1", "t2", "t3"),
    response_id = c("t1", "t2", "t3")
  )
  discovery_pairs <- data.frame(
    grna_target = c("t1", "t1", "t2", "t2"),
    response_id = c("response_4", "response_5", "response_4", "response_6")
  )

  scep <- import_data(
    grna_matrix = grna_matrix,
    response_matrix = response_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = "low"
  ) |>
    set_analysis_parameters(
      positive_control_pairs = positive_control_pairs,
      discovery_pairs = discovery_pairs,
      control_group = "nt_cells"
    ) |>
    assign_grnas(method = "thresholding", threshold = 40) |>
    # don't want to remove any cells or pairs for this one
    run_qc(
      response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0, 1),
      n_nonzero_trt_thresh = 0, n_nonzero_cntrl_thresh = 0
    ) |>
    run_calibration_check(calibration_group_size = 1) |>
    run_power_check() |>
    run_discovery_analysis()

  ## testing `plot_grna_count_distributions()` ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  plt <- plot_grna_count_distributions(scep, grnas_to_plot = c("id2", "id3", "nt2"))
  expect_equal(get_elements_from_plot_grid(plt, "title"), c("id2", "id3", "nt2"))
  expect_equal(get_elements_from_plot_grid(plt, "xlab-b"), c("gRNA count", "gRNA count", "gRNA count"))
  expect_equal(get_elements_from_plot_grid(plt, "ylab-l"), "Cell count (log scale)")

  ## testing `plot_assign_grnas()` ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pltlist <- plot_assign_grnas(scep, grnas_to_plot = c("nt1", "id1"), return_indiv_plots = TRUE)
  expect_equal(length(pltlist), 3)

  expect_true(!"title" %in% names(pltlist[[1]]$labels)) # no title for this plot
  expect_equal(pltlist[[1]]$labels$y, "gRNA count")
  expect_equal(pltlist[[1]]$labels$x, "Assignment")

  expect_equal(pltlist[[2]]$labels$title, "N cells per gRNA")
  expect_equal(pltlist[[2]]$labels$x, "N cells")
  expect_equal(pltlist[[2]]$labels$y, "gRNA")

  expect_equal(pltlist[[3]]$labels$title, "N gRNAs per cell")
  expect_equal(pltlist[[3]]$labels$x, "N gRNAs")
  expect_equal(pltlist[[3]]$labels$y, "Frequency")

  ## testing `plot_run_calibration_check()` ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pltlist <- plot_run_calibration_check(scep, return_indiv_plots = TRUE)
  expect_equal(length(pltlist), 4)

  expect_equal(pltlist[[1]]$labels$title, "QQ plot (bulk)")
  expect_equal(pltlist[[1]]$labels$x, "Expected null p-value")
  expect_equal(pltlist[[1]]$labels$y, "Observed p-value")

  expect_equal(pltlist[[2]]$labels$title, "QQ plot (tail)")
  expect_equal(pltlist[[2]]$labels$x, "Expected null p-value")
  expect_equal(pltlist[[2]]$labels$y, "Observed p-value")

  expect_equal(pltlist[[3]]$labels$title, "Log fold changes")
  expect_equal(pltlist[[3]]$labels$x, "Estimated log-2 fold change")
  expect_equal(pltlist[[3]]$labels$y, "Density")

  # this is the text-only pane
  expect_true(!"title" %in% names(pltlist[[4]]))
  expect_equal(pltlist[[4]]$labels$x, "x")
  expect_equal(pltlist[[4]]$labels$y, "y")

  # testing `plot_run_discovery_analysis()` ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pltlist <- plot_run_discovery_analysis(scep, return_indiv_plots = TRUE)
  expect_equal(length(pltlist), 4)

  expect_equal(pltlist[[1]]$labels$title, "QQ plot (bulk)")
  expect_equal(pltlist[[1]]$labels$x, "Expected null p-value")
  expect_equal(pltlist[[1]]$labels$y, "Observed p-value")

  expect_equal(pltlist[[2]]$labels$title, "QQ plot (tail)")
  expect_equal(pltlist[[2]]$labels$x, "Expected null p-value")
  expect_equal(pltlist[[2]]$labels$y, "Observed p-value")

  expect_equal(pltlist[[3]]$labels$title, "Discovery volcano plot")
  expect_equal(pltlist[[3]]$labels$x, "Log fold change")
  expect_equal(pltlist[[3]]$labels$y, "P-value")

  # this is the text-only pane
  expect_true(!"title" %in% names(pltlist[[4]]))
  expect_equal(pltlist[[4]]$labels$x, "x")
  expect_equal(pltlist[[4]]$labels$y, "y")

  # testing `plot_run_qc()` ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pltlist <- plot_run_qc(scep, return_indiv_plots = TRUE)
  expect_equal(length(pltlist), 2)

  expect_equal(pltlist[[1]]$labels$title, "Cellwise QC")
  expect_equal(pltlist[[1]]$labels$x, "Filter")
  expect_equal(pltlist[[1]]$labels$y, "Percent cells removed")

  expect_equal(pltlist[[2]]$labels$title, "Pairwise QC")
  expect_equal(pltlist[[2]]$labels$x, "N nonzero trt. cells")
  expect_equal(pltlist[[2]]$labels$y, "N nonzero cntrl. cells")

  # testing `plot_run_power_check()` ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  plt <- plot_run_power_check(scep)
  expect_equal(plt$labels$title, "Positive and negative control p-values")
  expect_equal(plt$labels$x, "Pair type")
  expect_equal(plt$labels$y, "p-value")

  # testing `plot_run_power_check()` ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  plt <- plot_response_grna_target_pair(scep, response_id = "response_4", grna_target = "t1")
  expect_equal(plt$labels$title, "Response: response_4\ngRNA target: t1")
  expect_equal(plt$labels$x, "Treatment status")
  expect_equal(plt$labels$y, "Normalized expression")
})
