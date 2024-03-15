test_that(".make_mock_target_and_chr_only", {
  inputs <- list(
    c(1, 1, 3),
    c(2, 5),
    c(20, 1),
    c(1, 1, 1),
    c(100, 10, 1),
    10,
    1
  )

  for (input in inputs) {
    out <- .make_mock_target_and_chr_only(num_guides_per_target = input)
    expect_equal(nrow(out), sum(input) * 2)
    expect_equal(length(unique(out$chr)), length(input) + 1)
    expect_equal(length(unique(out$grna_target)), length(input) * 2)
    expect_equal(strsplit(out$grna_target, "_") |> sapply(`[`, 2), out$chr)
    expect_equal(
      strsplit(out$grna_target, "_") |> sapply(`[`, 1),
      # in the first half each target is on its own chromosome so all are "t1"
      # in the second half all are on the same chromosome so the count increments
      rep(c("t1", paste0("t", 1:length(input))), c(sum(input), input))
    )
  }
})

test_that(".make_mock_cross_target_chr_with_distances", {
  target_and_chr_data_inputs <- list(
    c(1, 1, 3), c(2, 5), c(20, 1), 10
  )
  chr_distance_and_start_inputs <- list(
    list(distance = 0, start = 2),
    list(distance = 0:1, start = c(5, 1)),
    list(distance = 5, start = 0),
    list(distance = c(0, 8, 0, 1), start = c(1, 4, 2, 1))
  )

  for (target_and_chr_input in target_and_chr_data_inputs) {
    for (chr_distance_and_start_input in chr_distance_and_start_inputs) {
      target_and_chr <- .make_mock_target_and_chr_only(target_and_chr_input)

      out <- .make_mock_cross_target_chr_with_distances(
        target_and_chr_data = target_and_chr, chr_distances = chr_distance_and_start_input$distance,
        chr_starts = chr_distance_and_start_input$start
      )

      expect_equal(nrow(out), sum(target_and_chr_input) * 2 * length(chr_distance_and_start_input$distance))
      expect_equal(length(unique(out$chr)), (length(target_and_chr_input) + 1) * length(chr_distance_and_start_input$distance))
      expect_equal(length(unique(out$grna_target)), length(target_and_chr_input) * 2 * length(chr_distance_and_start_input$distance))
      # do target chr and actual chr match?
      expect_equal(
        strsplit(out$grna_target, "_") |> sapply(`[`, 2),
        strsplit(out$chr, "_") |> sapply(`[`, 1)
      )
      # do `end` and `start` agree with what's in `grna_target`?
      target_inds <- gsub(".*d", "", out$grna_target) |> as.numeric()
      expect_equal(out$start, chr_distance_and_start_input$start[target_inds])
      expect_equal(out$end - out$start, chr_distance_and_start_input$distance[target_inds])

      target_and_chr_target_starts <- rep(
        c("t1", paste0("t", 1:length(target_and_chr_input))),
        c(sum(target_and_chr_input), target_and_chr_input)
      )
      expect_equal(
        strsplit(out$grna_target, "_") |> sapply(`[`, 1),
        # in the first half of each `target_and_chr`, each target is on its own chromosome so all are "t1"
        # in the second half all are on the same chromosome so the count increments.
        # And then this is repeated once per each entry of chr_distance_and_start_input
        rep(target_and_chr_target_starts, times = length(chr_distance_and_start_input$distance))
      )
    }
  }
})

test_that("make_mock_grna_target_data", {
  target_and_chr_data_inputs <- list(
    c(1, 1, 3), c(2, 5), c(20, 1), 10
  )
  chr_distance_and_start_inputs <- list(
    list(distance = 0, start = 2),
    list(distance = 0:1, start = 5),
    list(distance = 5, start = 0),
    list(distance = c(0, 8, 0, 1), start = c(1, 4, 2, 1))
  )
  num_nts <- c(0, 5)

  for (target_and_chr_input in target_and_chr_data_inputs) {
    for (chr_distance_and_start_input in chr_distance_and_start_inputs) {
      for (num_nt in num_nts) {
        out <- make_mock_grna_target_data(
          num_guides_per_target = target_and_chr_input,
          chr_distances = chr_distance_and_start_input$distance,
          chr_starts = chr_distance_and_start_input$start,
          num_nt_guides = num_nt
        )
        expect_equal(length(unique(out$grna_id)), nrow(out))

        # `grna_target`, `chr`, `start` and `end` aren't modified
        # so this just tests `grna_id` and the appended NT guides
        out_pieces <- split(out, ifelse(out$grna_target == "non-targeting", "nt", "non_nt"))

        if (num_nt > 0) {
          expect_equal(nrow(out_pieces$nt), num_nt)
          expect_equal(sub("nt", "", out_pieces$nt$grna_id) |> as.numeric(), 1:num_nt)
        }
        expect_equal(
          sub("g.*_(t)", "\\1", out_pieces$non_nt$grna_id),
          out_pieces$non_nt$grna_target
        )
        with(out_pieces$non_nt, {
          for (target in unique(grna_target)) {
            id_start_nums <- grna_id[grna_target == target] |>
              strsplit("_") |>
              sapply(`[`, 1) |>
              sub(pattern = "g", replacement = "") |>
              as.numeric()
            expect_equal(id_start_nums, 1:sum(grna_target == target))
          }
        })
      }
    }
  }
})

test_that(".make_mock_patterned_matrix", {
  big_test <- 123

  for (num_rows in c(1, 5, 12, 13, 14, 25, 26, 27)) {
    for (num_cols in c(1, 5, 12, 13, 14, 25, 26, 27)) {
      out_cols <- .make_mock_patterned_matrix(num_rows, num_cols, patterns_at_col_level = TRUE, big = big_test)
      out_rows <- .make_mock_patterned_matrix(num_rows, num_cols, patterns_at_col_level = FALSE, big = big_test)

      expect_equal(dim(out_cols), c(num_rows, num_cols))
      expect_equal(dim(out_rows), c(num_rows, num_cols))
      expect_equal(out_cols[, 1], rep(0, num_rows))
      expect_equal(out_rows[1, ], rep(0, num_cols))

      ## col tests
      if (num_cols >= 2) {
        expect_equal(out_cols[, 2], rep(1, num_rows))
      }
      if (num_cols >= 3) {
        expect_equal(out_cols[, 3], rep(big_test, num_rows))
      }
      if (num_cols >= 4) {
        expect_equal(out_cols[, 4], rep(c(0, 1), c(num_rows - 1, 1)))
      }
      if (num_cols >= 5) {
        expect_equal(out_cols[, 5], rep(c(0, 1), c(1, num_rows - 1)))
      }
      if (num_cols >= 6) {
        expect_equal(out_cols[, 6], rep(c(big_test, 0), c(num_rows - 1, 1)))
      }
      if (num_cols >= 7) {
        expect_equal(out_cols[, 7], rep(c(big_test, 1), c(num_rows - 1, 1)))
      }
      if (num_cols >= 8) {
        expect_equal(out_cols[, 8], rep(c(0, big_test), c(num_rows - 1, 1)))
      }
      if (num_cols >= 9) {
        expect_equal(out_cols[, 9], rep(c(1, big_test), c(num_rows - 1, 1)))
      }
      if (num_cols >= 10) {
        expect_equal(out_cols[, 10], 0:(num_rows - 1))
      }
      if (num_cols >= 11) {
        expect_equal(out_cols[, 11], num_rows:1)
      }
      if (num_cols >= 12) {
        expect_equal(out_cols[, 12], seq(0, big_test, length = num_rows) |> round())
      }
      if (num_cols >= 13) {
        expect_equal(out_cols[, 13], seq(big_test, 1, length = num_rows) |> round())
      }
      if (num_cols > 13 && num_cols < 26) {
        # Prob(Pois(1) >= 20) < 2 * 10^{-19} so this won't be happening in these samples
        # even though >= 20 is a positive probability event (and with these seeds).
        # This is checking that `big` didn't sneak in there
        expect_true(all(out_cols[, 14:num_cols] < 20))
      }
      if (num_cols >= 26) {
        expect_equal(out_cols[, 1:13], out_cols[, 14:26])
      }

      ## row tests
      if (num_rows >= 2) {
        expect_equal(out_rows[2, ], rep(1, num_cols))
      }
      if (num_rows >= 3) {
        expect_equal(out_rows[3, ], rep(big_test, num_cols))
      }
      if (num_rows >= 4) {
        expect_equal(out_rows[4, ], rep(c(0, 1), c(num_cols - 1, 1)))
      }
      if (num_rows >= 5) {
        expect_equal(out_rows[5, ], rep(c(0, 1), c(1, num_cols - 1)))
      }
      if (num_rows >= 6) {
        expect_equal(out_rows[6, ], rep(c(big_test, 0), c(num_cols - 1, 1)))
      }
      if (num_rows >= 7) {
        expect_equal(out_rows[7, ], rep(c(big_test, 1), c(num_cols - 1, 1)))
      }
      if (num_rows >= 8) {
        expect_equal(out_rows[8, ], rep(c(0, big_test), c(num_cols - 1, 1)))
      }
      if (num_rows >= 9) {
        expect_equal(out_rows[9, ], rep(c(1, big_test), c(num_cols - 1, 1)))
      }
      if (num_rows >= 10) {
        expect_equal(out_rows[10, ], 0:(num_cols - 1))
      }
      if (num_rows >= 11) {
        expect_equal(out_rows[11, ], num_cols:1)
      }
      if (num_rows >= 12) {
        expect_equal(out_rows[12, ], seq(0, big_test, length = num_cols) |> round())
      }
      if (num_rows >= 13) {
        expect_equal(out_rows[13, ], seq(big_test, 1, length = num_cols) |> round())
      }
      if (num_rows > 13 && num_rows < 26) {
        # Prob(Pois(1) >= 20) < 2 * 10^{-19} so this won't be happening in these samples
        # even though >= 20 is a positive probability event (and with these seeds)
        expect_true(all(out_rows[14:num_rows, ] < 20))
      }
      if (num_rows >= 26) {
        expect_equal(out_rows[1:13, ], out_rows[14:26, ])
      }
    }
  }
})

test_that(".make_mock_matrix_pattern_list", {
  big <- 123
  dims_to_test <- list(c(1, 1), c(12, 25))
  for (dims in dims_to_test) {
    curr_nrow <- dims[1]
    curr_ncol <- dims[2]
    curr_rownames <- paste0("aaa_", 1:curr_nrow)

    ## testing just using "zero"
    results_zero <- .make_mock_matrix_pattern_list(
      patterns_to_make = "zero", pattern_nrow = curr_nrow, pattern_ncol = curr_ncol,
      big = big, pattern_rownames = curr_rownames
    )
    expect_equal(length(results_zero), 1)
    expect_equal(names(results_zero), "zero")
    expect_equal(results_zero[[1]], matrix(0, curr_nrow, curr_ncol, dimnames = list(curr_rownames)))

    ## testing using "one" and "big"
    results_one_and_big <- .make_mock_matrix_pattern_list(
      patterns_to_make = c("one", "big"), pattern_nrow = curr_nrow, pattern_ncol = curr_ncol,
      big = big, pattern_rownames = curr_rownames
    )
    expect_equal(length(results_one_and_big), 2)
    expect_equal(names(results_one_and_big), c("one", "big"))
    expect_equal(
      results_one_and_big[[1]],
      matrix(1, curr_nrow, curr_ncol, dimnames = list(curr_rownames))
    )
    expect_equal(
      results_one_and_big[[2]],
      matrix(big, curr_nrow, curr_ncol, dimnames = list(curr_rownames))
    )

    ## testing using all values
    set.seed(123) # for the calls to `.make_mock_patterned_matrix`
    results_all <- .make_mock_matrix_pattern_list(
      patterns_to_make = c("column", "row", "zero", "big", "one"),
      pattern_nrow = curr_nrow, pattern_ncol = curr_ncol,
      big = big, pattern_rownames = curr_rownames
    )
    expect_equal(length(results_all), 5)
    expect_equal(names(results_all), c("zero", "one", "big", "row", "column"))
    expect_equal(
      results_all[[1]],
      matrix(0, curr_nrow, curr_ncol, dimnames = list(curr_rownames))
    )
    expect_equal(
      results_all[[2]],
      matrix(1, curr_nrow, curr_ncol, dimnames = list(curr_rownames))
    )
    expect_equal(
      results_all[[3]],
      matrix(big, curr_nrow, curr_ncol, dimnames = list(curr_rownames))
    )
    set.seed(123)
    expect_equal(
      results_all[[4]],
      .make_mock_patterned_matrix(
        curr_nrow, curr_ncol,
        patterns_at_col_level = FALSE, big = big
      ) |> `rownames<-`(curr_rownames)
    )
    expect_equal(
      results_all[[5]],
      .make_mock_patterned_matrix(
        curr_nrow, curr_ncol,
        patterns_at_col_level = TRUE, big = big
      ) |> `rownames<-`(curr_rownames)
    )
  }
})


test_that("make_mock_grna_matrices", {
  grna_target_data_frame_no_nt <- make_mock_grna_target_data(
    num_guides_per_target = c(2, 3), chr_distances = 2,
    chr_starts = 1, num_nt_guides = 0
  )
  grna_target_data_frame_with_nt <- make_mock_grna_target_data(
    num_guides_per_target = 10, chr_distances = 2,
    chr_starts = 1, num_nt_guides = 5
  )
  num_cells <- 12
  big <- 123

  ##### 1. testing input processing of `non_nt_patterns` and `nt_patterns`
  expect_error(
    make_mock_grna_matrices(
      grna_target_data_frame = grna_target_data_frame_no_nt,
      num_cells = num_cells,
      non_nt_patterns = c("aaa", "zero") # this should cause an error
    ),
    regex = "`non_nt_patterns` must be a non-empty subset of"
  )
  expect_error(
    make_mock_grna_matrices(
      grna_target_data_frame = grna_target_data_frame_with_nt,
      num_cells = num_cells,
      non_nt_patterns = c("zero", "all")
      # the default value of `nt_patterns` should cause an error
    ),
    regex = "`grna_target_data_frame` contains non-targeting gRNA but `nt_patterns` has not been set"
  )
  expect_error(
    make_mock_grna_matrices(
      grna_target_data_frame = grna_target_data_frame_with_nt,
      num_cells = num_cells,
      non_nt_patterns = c("zero", "all"),
      nt_patterns = c("zero", "aaa") # this should cause an error
    ),
    regex = "`nt_patterns` must be a non-empty subset of"
  )

  ##### 2. return dimensions and names
  results_no_nt_one <- make_mock_grna_matrices(
    grna_target_data_frame = grna_target_data_frame_no_nt,
    num_cells = num_cells,
    non_nt_patterns = "one",
    big = big
  )
  set.seed(321)
  results_with_nt_single <- make_mock_grna_matrices(
    grna_target_data_frame = grna_target_data_frame_with_nt,
    num_cells = num_cells,
    non_nt_patterns = "column",
    nt_patterns = "row",
    big = big
  )

  set.seed(123)
  results_no_nt_several <- make_mock_grna_matrices(
    grna_target_data_frame = grna_target_data_frame_no_nt,
    num_cells = num_cells,
    non_nt_patterns = c("column", "row", "one"),
    big = big
  )
  results_with_nt_all <- make_mock_grna_matrices(
    grna_target_data_frame = grna_target_data_frame_with_nt,
    num_cells = num_cells,
    non_nt_patterns = "all",
    nt_patterns = "all",
    big = big
  )

  ## (a) checking `results_no_nt_one` which is a sparse matrix
  expect_true(is(results_no_nt_one, "TsparseMatrix"))
  expect_equal(dim(results_no_nt_one), c(nrow(grna_target_data_frame_no_nt), num_cells))
  expect_true(identical(rownames(results_no_nt_one), grna_target_data_frame_no_nt$grna_id))
  expect_equal(results_no_nt_one@x |> unique(), 1)

  ## (b) checking `results_with_nt_single` which is also a sparse matrix
  expect_true(is(results_with_nt_single, "TsparseMatrix"))
  expect_equal(dim(results_with_nt_single), c(nrow(grna_target_data_frame_with_nt), num_cells))
  expect_true(identical(rownames(results_with_nt_single), grna_target_data_frame_with_nt$grna_id))
  set.seed(321)
  expect_equal(
    results_with_nt_single[grna_target_data_frame_with_nt$grna_target != "non-targeting", ] |>
      as.matrix() |>
      `attr<-`("dimnames", NULL),
    .make_mock_patterned_matrix(sum(grna_target_data_frame_with_nt$grna_target != "non-targeting"), num_cells,
      patterns_at_col_level = TRUE, big = big
    )
  )
  set.seed(321)
  # this is just to get the seed right
  .make_mock_patterned_matrix(sum(grna_target_data_frame_with_nt$grna_target != "non-targeting"), num_cells,
    patterns_at_col_level = TRUE, big = big
  )
  expect_equal(
    results_with_nt_single[grna_target_data_frame_with_nt$grna_target == "non-targeting", ] |>
      as.matrix() |>
      `attr<-`("dimnames", NULL),
    .make_mock_patterned_matrix(sum(grna_target_data_frame_with_nt$grna_target == "non-targeting"), num_cells,
      patterns_at_col_level = FALSE, big = big
    )
  )

  ## (c) checking `results_no_nt_all` which is a list of 5 sparse matrices
  expect_true(is(results_no_nt_several, "list"))
  expect_true(all(sapply(results_no_nt_several, is, "TsparseMatrix")))
  expect_equal(length(results_no_nt_several), 3)
  expect_equal(names(results_no_nt_several), paste0("non_nt_", c("one", "row", "column")))

  expect_true(all(sapply(results_no_nt_several, function(mat) dim(mat) == c(nrow(grna_target_data_frame_no_nt), num_cells))))
  expect_true(all(sapply(results_no_nt_several, function(mat) identical(rownames(mat), grna_target_data_frame_no_nt$grna_id))))
  set.seed(123)
  expect_equal(
    results_no_nt_several$non_nt_row |> as.matrix() |> `attr<-`("dimnames", NULL),
    .make_mock_patterned_matrix(nrow(grna_target_data_frame_no_nt), num_cells,
      patterns_at_col_level = FALSE, big = big
    )
  )

  ## (d) checking `results_with_nt_all` which is a list of 25 sparse matrices
  expect_true(is(results_with_nt_all, "list"))
  expect_true(all(sapply(results_with_nt_all, is, "TsparseMatrix")))
  expect_equal(length(results_with_nt_all), 25)

  pattern_names <- c("zero", "one", "big", "row", "column")
  all_names <- outer(paste0("non_nt_", pattern_names, "_"), paste0("nt_", pattern_names), paste0) |>
    t() |>
    as.character()

  expect_equal(names(results_with_nt_all), all_names)

  expect_true(all(sapply(results_with_nt_all, function(mat) dim(mat) == c(nrow(grna_target_data_frame_with_nt), num_cells))))
  expect_true(all(sapply(results_with_nt_all, function(mat) identical(rownames(mat), grna_target_data_frame_with_nt$grna_id))))

  ##### 3. testing `unlist_if_single_pattern` and `return_as_sparse`
  return_no_unlist <- make_mock_grna_matrices(
    grna_target_data_frame = grna_target_data_frame_no_nt,
    num_cells = num_cells,
    non_nt_patterns = "one",
    nt_patterns = "big",
    unlist_if_single_pattern = FALSE,
    return_as_sparse = FALSE
  )
  expect_true(is(return_no_unlist, "list"))
  expect_equal(length(return_no_unlist), 1)
  expect_false(is(return_no_unlist[[1]], "TsparseMatrix"))
  expect_true(all(return_no_unlist[[1]] %in% c(1, big)))
})


test_that("make_mock_response_matrices", {
  num_responses <- 18
  num_cells <- 12
  big <- 123

  ## testing input processing of `patterns`
  expect_error(
    make_mock_response_matrices(
      num_responses = num_responses,
      num_cells = num_cells,
      patterns = c("aaa", "zero") # this should cause an error
    ),
    regex = "`patterns` must be a non-empty subset of"
  )

  ## testing single return
  set.seed(191)
  results_single <- make_mock_response_matrices(
    num_responses = num_responses,
    num_cells = num_cells,
    patterns = "row",
    return_as_sparse = FALSE,
    big = big
  )

  expect_true(is(results_single, "matrix"))
  set.seed(191)
  expect_equal(
    results_single |> as.matrix() |> `attr<-`("dimnames", NULL),
    .make_mock_patterned_matrix(num_responses, num_cells, patterns_at_col_level = FALSE, big = big)
  )

  ## testing list return
  set.seed(919)
  results_all <- make_mock_response_matrices(
    num_responses = num_responses,
    num_cells = num_cells,
    patterns = "all",
    big = big
  )

  expect_true(is(results_all, "list"))
  expect_equal(length(results_all), 5)
  expect_equal(names(results_all), c("zero", "one", "big", "row", "column"))
  set.seed(919)
  expect_equal(
    results_all$row |> as.matrix() |> `attr<-`("dimnames", NULL),
    .make_mock_patterned_matrix(num_responses, num_cells, patterns_at_col_level = FALSE, big = big)
  )
})

test_that("make_mock_extra_covariates_data_frames", {
  # this sequence is chosen so that we have values on either side of the various
  # values that divisibility checks are used for
  num_cell_values <- c(20, 21, 22, 23, 24, 25, 26, 49, 50, 51)

  all_names <- c(
    "constant", "almost_constant_two_levels_one_val",
    "almost_constant_two_levels_two_vals", "almost_constant_many_levels_one_val",
    "almost_constant_many_levels_two_vals", "many_levels", "missing_level",
    "numeric", "count", "factor_one_value_level", "factor_many_levels", "many_columns"
  )

  ## testing results when `patterns = "all"` so that all computations are verified
  for (num_cells in num_cell_values) {
    set.seed(num_cells)
    results_all <- make_mock_extra_covariates_data_frames(num_cells, patterns = "all")

    expect_length(results_all, 12)
    expect_equal(names(results_all), all_names)
    # first 7 elements just test `batch`, then next 4 elements test adding a second column
    # and last one tests adding 3 columns to `batch`
    expect_equal(sapply(results_all, ncol) |> as.numeric(), c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4))
    expect_equal(sapply(results_all, nrow) |> as.numeric(), rep(num_cells, 12))

    expect_equal(results_all$constant$batch, rep("b1", num_cells) |> factor(levels = "b1"))
    expect_equal(table(results_all$almost_constant_two_levels_one_val$batch) |> as.numeric(), c(1, num_cells - 1))
    expect_equal(table(results_all$almost_constant_two_levels_two_vals$batch) |> as.numeric(), c(2, num_cells - 2))
    expect_equal(table(results_all$almost_constant_many_levels_one_val$batch) |> as.numeric(), rep(c(1, num_cells - 9), c(9, 1)))
    expect_equal(table(results_all$almost_constant_many_levels_two_vals$batch) |> as.numeric(), rep(c(2, num_cells - 18), c(9, 1)))
    expect_equal(table(results_all$many_levels$batch) |> as.numeric(), rep(c(num_cells %/% 5, num_cells %/% 5 + num_cells %% 5), c(4, 1)))
    expect_equal(table(results_all$missing_level$batch) |> as.numeric(), c(num_cells %/% 2, num_cells %/% 2 + num_cells %% 2, 0))

    expect_equal(table(results_all$factor_one_value_level$batch) |> as.numeric(), c(num_cells %/% 3, num_cells %/% 3, num_cells %/% 3 + num_cells %% 3))
    expect_equal(table(results_all$factor_many_levels$factor) |> as.numeric(), num_cells %/% 5 + (num_cells %% 5) * rep(0:1, c(4, 1)))
    expect_equal(table(results_all$many_columns$factor) |> as.numeric(), num_cells %/% 5 + (num_cells %% 5) * rep(0:1, c(4, 1)))
  }

  ## testing single return using last value of `num_cells` from the above loop
  set.seed(num_cells)
  results_single_df <- make_mock_extra_covariates_data_frames(num_cells, patterns = "many_columns")
  expect_true(is(results_single_df, "data.frame"))
  # same as what we get when more patterns are used?
  set.seed(num_cells)
  results_several_df <- make_mock_extra_covariates_data_frames(num_cells, patterns = c("many_columns", "constant", "almost_constant_many_levels_two_vals"))
  expect_equal(results_single_df, results_several_df[[3]])
})
