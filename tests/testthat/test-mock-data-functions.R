test_that(".make_mock_target_and_chr_only", {
  inputs <- list(
    c(1, 1, 3),
    c(2, 5),
    c(20, 1),
    c(1,1,1),
    c(100,10,1),
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
    list(distance = 0:1, start = c(5,1)),
    list(distance = 5, start = 0),
    list(distance = c(0,8,0,1), start = c(1,4,2,1))
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
      expect_equal(strsplit(out$grna_target, "_") |> sapply(`[`, 2),
                   strsplit(out$chr, "_") |> sapply(`[`, 1))
      # do `end` and `start` agree with what's in `grna_target`?
      target_inds <- gsub(".*d", "", out$grna_target) |> as.numeric()
      expect_equal(out$start, chr_distance_and_start_input$start[target_inds])
      expect_equal(out$end - out$start, chr_distance_and_start_input$distance[target_inds])

      target_and_chr_target_starts <- rep(c("t1", paste0("t", 1:length(target_and_chr_input))),
                                          c(sum(target_and_chr_input), target_and_chr_input))
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
    list(distance = c(0,8,0,1), start = c(1,4,2,1))
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
        expect_equal(sub("g.*_(t)", "\\1", out_pieces$non_nt$grna_id),
                     out_pieces$non_nt$grna_target)
        with(out_pieces$non_nt, {
          for (target in unique(grna_target)) {
            id_start_nums <- grna_id[grna_target == target] |>
              strsplit("_") |> sapply(`[`, 1) |> sub(pattern = "g", replacement = "") |>
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
      out_cols <- .make_mock_patterned_matrix(num_rows, num_cols, patterns_at_col_level = TRUE, big = big_test, seed = num_rows + num_cols)
      out_rows <- .make_mock_patterned_matrix(num_rows, num_cols, patterns_at_col_level = FALSE, big = big_test, seed = num_rows + num_cols)

      expect_equal(dim(out_cols), c(num_rows, num_cols))
      expect_equal(dim(out_rows), c(num_rows, num_cols))
      expect_equal(out_cols[,1], rep(0, num_rows))
      expect_equal(out_rows[1,], rep(0, num_cols))

      ## col tests
      if (num_cols >= 2) {
        expect_equal(out_cols[,2], rep(1, num_rows))
      }
      if (num_cols >= 3) {
        expect_equal(out_cols[,3], rep(big_test, num_rows))
      }
      if (num_cols >= 4) {
        expect_equal(out_cols[,4], rep(c(0, 1), c(num_rows- 1, 1)))
      }
      if (num_cols >= 5) {
        expect_equal(out_cols[,5], rep(c(0, 1), c(1, num_rows - 1)))
      }
      if (num_cols >= 6) {
        expect_equal(out_cols[,6], rep(c(big_test, 0), c(num_rows - 1, 1)))
      }
      if (num_cols >= 7) {
        expect_equal(out_cols[,7], rep(c(big_test, 1), c(num_rows - 1, 1)))
      }
      if (num_cols >= 8) {
        expect_equal(out_cols[,8], rep(c(0, big_test), c(num_rows - 1, 1)))
      }
      if (num_cols >= 9) {
        expect_equal(out_cols[,9], rep(c(1, big_test), c(num_rows - 1, 1)))
      }
      if (num_cols >= 10) {
        expect_equal(out_cols[,10], 0:(num_rows-1))
      }
      if (num_cols >= 11) {
        expect_equal(out_cols[,11], num_rows:1)
      }
      if (num_cols >= 12) {
        expect_equal(out_cols[,12], seq(0, big_test, length = num_rows) |> round())
      }
      if (num_cols >= 13) {
        expect_equal(out_cols[,13], seq(big_test, 1, length = num_rows) |> round())
      }
      if (num_cols > 13 && num_cols < 26) {
        # Prob(Pois(1) >= 20) < 2 * 10^{-19} so this won't be happening in these samples
        # even though >= 20 is a positive probability event (and with these seeds).
        # This is checking that `big` didn't sneak in there
        expect_true(all(out_cols[,14:num_cols] < 20))
      }
      if (num_cols >= 26) {
        expect_equal(out_cols[,1:13], out_cols[,14:26])
      }

      ## row tests
      if (num_rows >= 2) {
        expect_equal(out_rows[2,], rep(1, num_cols))
      }
      if (num_rows >= 3) {
        expect_equal(out_rows[3,], rep(big_test, num_cols))
      }
      if (num_rows >= 4) {
        expect_equal(out_rows[4,], rep(c(0, 1), c(num_cols- 1, 1)))
      }
      if (num_rows >= 5) {
        expect_equal(out_rows[5,], rep(c(0, 1), c(1, num_cols - 1)))
      }
      if (num_rows >= 6) {
        expect_equal(out_rows[6,], rep(c(big_test, 0), c(num_cols - 1, 1)))
      }
      if (num_rows >= 7) {
        expect_equal(out_rows[7,], rep(c(big_test, 1), c(num_cols - 1, 1)))
      }
      if (num_rows >= 8) {
        expect_equal(out_rows[8,], rep(c(0, big_test), c(num_cols - 1, 1)))
      }
      if (num_rows >= 9) {
        expect_equal(out_rows[9,], rep(c(1, big_test), c(num_cols - 1, 1)))
      }
      if (num_rows >= 10) {
        expect_equal(out_rows[10,], 0:(num_cols-1))
      }
      if (num_rows >= 11) {
        expect_equal(out_rows[11,], num_cols:1)
      }
      if (num_rows >= 12) {
        expect_equal(out_rows[12,], seq(0, big_test, length = num_cols) |> round())
      }
      if (num_rows >= 13) {
        expect_equal(out_rows[13,], seq(big_test, 1, length = num_cols) |> round())
      }
      if (num_rows > 13 && num_rows < 26) {
        # Prob(Pois(1) >= 20) < 2 * 10^{-19} so this won't be happening in these samples
        # even though >= 20 is a positive probability event (and with these seeds)
        expect_true(all(out_rows[14:num_rows,] < 20))
      }
      if (num_rows >= 26) {
        expect_equal(out_rows[1:13,], out_rows[14:26,])
      }
    }
  }
})

test_that("make_mock_grna_matrix_list", {
  grna_target_data_frame_list <- list(
    make_mock_grna_target_data(
      num_guides_per_target = 10, chr_distances = 2,
      chr_starts = 1, num_nt_guides = 0
    ),
    make_mock_grna_target_data(
      num_guides_per_target = 10, chr_distances = 2,
      chr_starts = 1, num_nt_guides = 0
    ) |>
      dplyr::filter(grna_target == "t1_c1_d1"),  # taking just 1 target and chr
    make_mock_grna_target_data(
      num_guides_per_target = 10, chr_distances = 2,
      chr_starts = 1, num_nt_guides = 5
    ),
    make_mock_grna_target_data(
      num_guides_per_target = c(1,5,10), chr_distances = 2,
      chr_starts = 1, num_nt_guides = 0
    ),
    make_mock_grna_target_data(
      num_guides_per_target = c(5, 5, 5), chr_distances = 2,
      chr_starts = 1, num_nt_guides = 4
    ) |> dplyr::filter(chr == "c4_d1"), # taking one bigger chr
    make_mock_grna_target_data(
      num_guides_per_target = c(1,5,10), chr_distances = 2,
      chr_starts = 1, num_nt_guides = 4
    )
  )

  num_cell_values <- c(1, 5, 18)
  big_test <- 321

  for (grna_target_data_frame in grna_target_data_frame_list) {
    for (num_cells in num_cell_values) {
      out <- make_mock_grna_matrix_list(grna_target_data_frame, num_cells, big = big_test,
                                        seed = nrow(grna_target_data_frame) + num_cells, return_as_sparse = FALSE)

      expect_true(all(sapply(out, nrow) == nrow(grna_target_data_frame)))
      expect_true(all(sapply(out, ncol) == num_cells))
      expect_true(all(sapply(out, function(df) all(rownames(df) == grna_target_data_frame$grna_id))))

      num_nt <- sum(grna_target_data_frame$grna_target == "non-targeting")
      expect_length(out, ifelse(num_nt > 0, 25, 5))

      expect_true(all(out[[1]] == 0))

      if (num_nt == 0) {
        # then out 2:5 are the predicted patterns and that's all there is
        expect_true(all(out[[2]] == 1))
        expect_true(all(out[[3]] == big_test))
        expect_true(all(out[[4]] == .make_mock_patterned_matrix(nrow(grna_target_data_frame), num_cells, FALSE, big = big_test,
                                                                seed = nrow(grna_target_data_frame) + num_cells)))
        expect_true(all(out[[5]] == .make_mock_patterned_matrix(nrow(grna_target_data_frame), num_cells, TRUE, big = big_test,
                                                                seed = nrow(grna_target_data_frame) + num_cells)))
      } else {
        # the non-NT patterns are tested above so this section will just test the NT patterns
        # NT patterns are looped thru first so we just need out[2:5] again
        nt_submatrices <- lapply(out[2:5], function(m) m[grna_target_data_frame$grna_target == "non-targeting", ,drop = FALSE])
        expect_true(all(nt_submatrices[[1]] == 1))
        expect_true(all(nt_submatrices[[2]] == big_test))
        expect_true(all(nt_submatrices[[3]][1,] == seq(big_test, 1, length = num_cells) |> round()))
        expect_true(all(nt_submatrices[[4]] == .make_mock_patterned_matrix(num_nt, num_cells, TRUE, big = big_test,
                                                                           seed = nrow(grna_target_data_frame) + num_cells)))
      }
    }
  }
})

test_that("make_mock_response_matrix_list ", {

  num_response_values <- c(1, 4, 13, 14, 25, 26, 27)
  num_cell_values <- c(1, 5, 13, 14, 25, 26, 27)
  big_test <- 345

  for (num_responses in num_response_values) {
    for (num_cells in num_cell_values) {
      out <- make_mock_response_matrix_list(num_responses, num_cells, big = big_test,
                                        seed = num_response_values + num_cells, return_as_sparse = FALSE)

      expect_true(all(sapply(out, nrow) == num_responses))
      expect_true(all(sapply(out, ncol) == num_cells))
      expect_true(all(sapply(out, function(df) rownames(df) == paste0("response_", 1:num_responses))))
      # expect_true(all(sapply(out, function(df) colnames(df) == paste0("cell_", 1:num_cells))))  # not currently using cell names

      expect_true(all(out[[1]] == 0))
      expect_true(all(out[[2]] == 1))
      expect_true(all(out[[3]] == big_test))
      expect_true(all(out[[4]] == .make_mock_patterned_matrix(num_responses, num_cells, FALSE, big = big_test,
                                                              seed = num_response_values + num_cells)))
      expect_true(all(out[[5]] == .make_mock_patterned_matrix(num_responses, num_cells, TRUE, big = big_test,
                                                              seed = num_response_values + num_cells)))
    }
  }
})

test_that("make_mock_extra_covariates_list", {

  num_cell_values <- c(20, 21, 22, 23, 24, 25, 26, 49, 50, 51)

  for (num_cells in num_cell_values) {
    out <- make_mock_extra_covariates_list(num_cells)

    expect_length(out, 12)
    # first 7 elements just test `batch`, then next 4 elements test adding a second column
    # and last one tests adding 3 columns to `batch`
    expect_equal(sapply(out, ncol) |> as.numeric(), c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4))
    expect_equal(sapply(out, nrow) |> as.numeric(), rep(num_cells, 12))

    expect_equal(out$constant$batch, rep("b1", num_cells) |> factor(levels = "b1"))
    expect_equal(table(out$almost_constant_two_levels_one_val$batch) |> as.numeric(), c(1, num_cells - 1))
    expect_equal(table(out$almost_constant_two_levels_two_vals$batch) |> as.numeric(), c(2, num_cells - 2))
    expect_equal(table(out$almost_constant_many_levels_one_val$batch) |> as.numeric(), rep(c(1, num_cells - 9), c(9, 1)))
    expect_equal(table(out$almost_constant_many_levels_two_vals$batch) |> as.numeric(), rep(c(2, num_cells - 18), c(9, 1)))
    expect_equal(table(out$many_levels$batch) |> as.numeric(), rep(c(num_cells %/% 10, num_cells %/% 10 + num_cells %% 10), c(9, 1)))
    expect_equal(table(out$missing_level$batch) |> as.numeric(), c(num_cells %/% 2, num_cells %/% 2 + num_cells %% 2, 0))

    expect_equal(table(out$factor_one_value_level$batch) |> as.numeric(), c(num_cells %/% 3, num_cells %/% 3, num_cells %/% 3 + num_cells %% 3))
    expect_equal(table(out$factor_many_levels$factor) |> as.numeric(), num_cells %/% 5 + (num_cells %% 5) * rep(0:1, c(4,1)))
    expect_equal(table(out$many_columns$factor) |> as.numeric(), num_cells %/% 5 + (num_cells %% 5) * rep(0:1, c(4,1)))
  }
})
