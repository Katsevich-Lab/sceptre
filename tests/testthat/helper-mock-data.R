######################################################################################################
# This section has old functions that are being replaced. Do not use.

#' #' Makes a mock UMI count matrix for testing.
#' #'
#' #' This can replicate either a response matrix or a gRNA expression matrix.
#' #'
#' #' The specific patterns that appear:
#' #' * constant
#' #' * all but one value constant
#' #' * linear sequences
#' #' * random
#' #'
#' #' There are 13 deterministic patterns tested, so with at least one random pattern we
#' #' require at least 14 in the dimension in which the patterns go; in other words, if
#' #' \code{patterns_at_col_level = TRUE} then \code{num_cells >= 14} is required (since
#' #' the patterns are at the column / cell level), and if instead \code{patterns_at_col_level = FALSE}
#' #' then \code{num_genes >= 14} is required.
#' #'
#' #'
#' #' @param num_genes, @param num_cells : the number of rows and columns in the resulting data, respectively
#' #' @param rowname_base : the rownames are of the form `rowname_base` + "_" + (1:num_genes)
#' #' @param patterns_at_col_level : if TRUE the columns contain the above patterns, else the rows do
#' #' @param big : a constant representing a large data value
#' #' @param seed : for the random part of the data
#' #'
#' #' @return a sparse matrix of class \code{dgTMatrix} mocking a response matrix.
#' #'
#' #'
#' #' @examples
#' #' # mock response UMI matrix with 5 rows and 20 columns
#' #' # the structured patterns are for each cell, so e.g. cell 1 has UMI counts of all 0
#' #' make_mock_count_matrix(num_genes = 5, num_cells = 20, rowname_base = "gene")
#' #'
#' #' # mock response UMI matrix with 20 rows and 5 columns
#' #' # the structured patterns are now for each gene, so e.g. gene 1 has UMI counts of all 0
#' #' make_mock_count_matrix(num_genes = 20, num_cells = 5, rowname_base = "gene", patterns_at_col_level = FALSE)
#' #'
#' #' # we can also mock gRNA count matrices
#' #' make_mock_count_matrix(num_genes = 7, num_cells = 20, rowname_base = "grna")
#' make_mock_count_matrix <- function(num_genes, num_cells, rowname_base, patterns_at_col_level = TRUE, big = 1000, seed = NULL) {
#'   if(is.null(seed)) {
#'     seed = 101000
#'   }
#'   set.seed(seed)
#'   dims <- if(patterns_at_col_level) c(num_genes, num_cells) else c(num_cells, num_genes)
#'   if(dims[2] < 13) {
#'     stop("There are 13 required patterns and the provided dimensions are too small.")
#'   }
#'   mock_response_matrix <- cbind(
#'     0, 1, big, # constant columns
#'     rep(c(0,1), c(dims[1]-1, 1)),  rep(c(0, 1), c(1, dims[1]-1)), # all but one element constant
#'     rep(c(big, 0), c(dims[1]-1, 1)) ,rep(c(big, 1), c(dims[1]-1, 1)),
#'     rep(c(0, big), c(dims[1]-1, 1)), rep(c(1, big), c(dims[1]-1, 1)),
#'     0:(dims[1] - 1), dims[1]:1, # sequence columns
#'     (0:(dims[1] - 1)) * big, (dims[1]:1) * big
#'   ) |>
#'     `colnames<-`(NULL) # the use of `big` adds a single column name
#'   # if we have more columns than the deterministic patterns, then add in random features
#'   if(dims[2] >= 14) {
#'     mock_response_matrix <- cbind(
#'       mock_response_matrix,
#'       matrix(
#'         rpois(dims[1] * (dims[2] - ncol(mock_response_matrix)), 1),
#'         nrow = dims[1]
#'       )
#'     )
#'   }
#'   if(!patterns_at_col_level) {
#'     mock_response_matrix <- t(mock_response_matrix)
#'   }
#'   rownames(mock_response_matrix) <- paste0(rowname_base, "_", 1:nrow(mock_response_matrix))  # rows are always genes
#'   return(as(mock_response_matrix, "TsparseMatrix"))
#' }

# `sample` has unwanted behavior here where e.g. `sample(2:2, 1)` can return 1
# even though only returning 2 is desired.
sample_as_vec <- function(x) {
  if(length(x) == 1)
    return(x)
  else {
    return(sample(x, 1))
  }
}

#' Helper function to make a factor with random numbers of entries per level.
#'
#' There are guaranteed to be at least two entries per level.
#'
#' @param num_factor_levels : int, the number of distinct values for this factor to have
#' @param num_entries : int, the length of the resulting vector
#' @param level_name_base : string, all values start with this + "_".
#' @param return_factor : bool, if \code{FALSE} the result is a character vector while if \code{TRUE}
#' it is actually a factor.
#' @param shuffle : if TRUE then the values of the vector are shuffled
#'
#' @return
#'
#'
#' @examples
make_random_factor <- function(num_factor_levels, num_entries, level_name_base = "lev", return_factor = TRUE, shuffle = FALSE) {
  lev_counts <- numeric(num_factor_levels)
  lev_counts[1] <- sample_as_vec(2:(num_entries - 2 * (num_factor_levels - 1)))
  if(num_factor_levels > 2) { # middle levels get random counts
    for(j in 2:(num_factor_levels-1)) {
      # we need to save 2 * (# levels remaining) entries for the future levels, and
      # we've allocated `sum(lev_counts)` so far
      lev_counts[j] <- sample_as_vec(2:(num_entries - 2 * (num_factor_levels - j) - sum(lev_counts)))
    }
  }
  lev_counts[num_factor_levels] <- num_entries - sum(lev_counts) # needs to sum to `num_entries`
  lev_names <- paste0(level_name_base, "_", 1:num_factor_levels)
  values <- rep(lev_names, lev_counts)
  if(shuffle) {
    values <- sample(values, num_entries, FALSE)
  }
  if(return_factor) {
    values <- factor(values, levels = lev_names)
  }
  return(values)
}

#' Mock cell-level covariates
#'
#' @param rep_level_counts : this is for mocking a replicate or batch feature. \code{length(rep_level_counts)}
#' gives the number of replicates, and each entry is the number of cells within that replicate.
#' \code{sum(rep_level_counts)} is the total number of cells. For example,
#' \code{rep_level_counts = c(2,4,6)} will result in data for 12 cells total, with
#' replicate sizes of 2, 4, and 6 respectively.
#' @param add_numeric : integter, if positive then that many columns of
#' iid N(0,1) data are added, mocking continuous features.
#' @param add_count integer, if positive then that many columns of
#' iid Pois(1) data are added, mocking count features.
#' @param add_factor integer, if positive then that many random factor features
#' are added. For each random factor, the number of levels is sampled uniformly from
#' \code{2:max_num_factor_levels}, and then the number of cells per level are filled
#' in randomly, with the constraint that no level has fewer than 2 cells.
#'
#' @param max_num_factor_levels : integer, the maximum number of levels for each random
#' factor. Ignored if \code{add_factor = 0}.
#' @param seed
#'
#' @return a data.frame with the specified covariates, with each row corresponding to a cell.
#'
#'
#' @examples
#' # just doing replicate for 20 cells
#' make_mock_extra_covariates(rep_level_counts = c(4, 13, 3))
#' # simulating features of all types. With this seed, two of the factors ended up with
#' # just 2 levels, while the third ended up with 4 levels.
#' make_mock_extra_covariates(rep_level_counts = c(4, 13, 3), add_numeric = 2, add_count = 1,
#' add_factor = 3, max_num_factor_levels = 4, seed = 111)
make_mock_extra_covariates <- function(rep_level_counts, add_numeric = 0, add_count = 0, add_factor = 0,
                                       max_num_factor_levels = 6, seed = NULL) {
  if(is.null(seed)) {
    seed = 101001
  }
  set.seed(seed)

  num_cells <- sum(rep_level_counts)
  if(add_factor > 0 && max_num_factor_levels * 2 > num_cells) {
    stop("`max_num_factor_levels` must be decreased so that at least 2 cells per level are possible.")
  }
  if(max_num_factor_levels < 2) {
    stop("`max_num_factor_levels` must be at least 2.")
  }
  covariates <- data.frame(
    replicate = rep(paste0("rep_", 1:length(rep_level_counts)), rep_level_counts) |>
      factor(levels = paste0("rep_", 1:length(rep_level_counts)))
  )
  if(add_numeric > 0) {
    for(i in 1:add_numeric) {
      covariates[[paste0("numeric_", i)]] <- rnorm(num_cells)
    }
  }
  if(add_count > 0) {
    for(i in 1:add_count) {
      covariates[[paste0("count_", i)]] <- rpois(num_cells, 1)
    }
  }
  if(add_factor > 0) {
    for(i in 1:add_factor) {
      num_levs <- sample(2:max_num_factor_levels, 1)
      covariates[[paste0("factor_", i)]] <- make_random_factor(
        num_factor_levels = num_levs, num_entries = num_cells,
        level_name_base = paste0("factor_", i),
        return_factor = TRUE, shuffle = TRUE
      )
    }
  }
  return(covariates)
}

#############
## Example ## [UNDER CONSTRUCTION]
#############

# num_cells <- 15
# num_genes <- 18
# num_grnas <- 5
# num_distinct_targets <- 4
# rep_levels <- c(4, 8, 3); stopifnot(sum(rep_levels) == num_cells)
# (mock_response <- make_mock_count_matrix(num_genes = num_genes, num_cells = num_cells, rowname_base = "gene", patterns_at_col_level = FALSE))
# (mock_grna <- make_mock_count_matrix(num_genes = num_grnas, num_cells = num_cells, rowname_base = "grna", patterns_at_col_level = TRUE))
# (mock_covariates <- make_mock_extra_covariates(rep_level_counts = rep_levels))
# ## if we want a fancier one
# # make_mock_extra_covariates(rep_level_counts = rep_levels, add_numeric = 2, add_count = 1, add_factor = 3, max_num_factor_levels = 4)
# (mock_grna_target_data_frame <- make_mock_grna_target_data(rownames(mock_grna), num_distinct_targets))







######################################################################################################
##                          Functions for making mock data for unit testing                         ##
######################################################################################################


# This function takes in a vector like c(1,5,3) and returns a data.frame with
# columns "grna_target" and "chr" with rows representing different genomic targets
# and what chromosome they are on.
# The input determines how many different guides there will be for each target, so
# in the example `num_guides_per_target = c(1, 5, 3)` we will have one target with one guide,
# one target with 5 guides, and one target with 3 guides. These will all be on different chromosomes.
# This is then repeated except all are on the same chromosome, so this example will have 4 chromosomes, with
# the final chromosome having length(c(1, 5, 3)) targets with 1+5+3 guides total for it.
.make_mock_target_and_chr_only <- function(num_guides_per_target) {
  targets <- rep(paste0("t", 1:length(num_guides_per_target)), num_guides_per_target)
  chrs <- rep(paste0("c", 1:(length(num_guides_per_target) + 1)), c(num_guides_per_target, sum(num_guides_per_target)))
  grna_target_data_simple <- data.frame(
    grna_target = c(
      paste0("t1_c", rep(1:length(num_guides_per_target), num_guides_per_target)),
      paste0(targets, "_c", length(num_guides_per_target) + 1)
    ),
    chr = chrs
  )
  return(grna_target_data_simple)
}

# This function takes the output of .make_mock_target_and_chr_only and adds in columns "start" and "end"
# representing where the target starts and stops on the given chromosome.
# This function is not meant to be directly called. It is assumed `length(chr_distances) == length(chr_starts)`.
# The result is a data.frame with columns "grna_target", "chr", "start", and "end", consisting of
# `target_and_chr_data` copied once per entry of `chr_distances` with each copy having the corresponding
# values of `chr_distances` and `chr_starts`.
# The input's `grna_target` and `chr` are modified to indicate which distance they correspond to.
# Note: in the semantic naming of targets and chromosomes, the index of the distance is used rather
# than the value of the distance, because there could be repeated distances.
.make_mock_cross_target_chr_with_distances <- function(target_and_chr_data, chr_distances, chr_starts) {
  target_chr_and_dist_data <- lapply(1:length(chr_distances), function(i) {
    target_and_chr_data |>
      dplyr::transmute(
        grna_target = paste0(grna_target, "_d", i),
        chr = paste0(chr, "_d", i),
        start = chr_starts[i], end = chr_starts[i] + chr_distances[i]
      )
  }) |>
    do.call(what = rbind)
  return(target_chr_and_dist_data)
}

#' Make mock grna_target_data_frame's
#'
#' This function makes a mock grna_target_data_frame in the spirit of
#' \code{data(grna_target_data_frame_highmoi)} and \code{data(grna_target_data_frame_lowmoi)}.
#' It has columns \code{grna_id}, \code{grna_target}, \code{chr}, \code{start}, and \code{end}.
#'
#' It builds the mock data in three steps.
#'
#' The first step takes in \code{num_guides_per_target} and makes a data.frame with one grna target
#' per entry of \code{num_guides_per_target}, repeated according to the values in that vector. Each
#' of these grna targets is on its own chromosome. This data.frame is then repeated except all targets
#' are now on the same chromosome.
#'
#' The second step adds in \code{start} and \code{end}. For every unique value of \code{chr_distances}
#' (and the corresponding value of \code{chr_starts}), the previous step is duplicated and then this copy's
#' \code{start} and \code{end} are filled in with the corresponding values of \code{chr_distances} and
#' \code{chr_starts}.
#'
#' Finally, the \code{grna_id} column is added in and optionally NT guides are appended.
#'
#' @param num_guides_per_target : a vector (length >=1) of positive ints giving the number of distinct grnas per target,
#' before the doubling in the chromosome addition step and the multiplication in the \code{start} and \code{end}
#' steps as described above.
#' @param chr_distances : a vector (length >=1) of non-negative ints giving the length on the chromosome of the various
#' guides.
#' @param chr_starts : a vector of positive ints either of length \code{chr_distances} or length 1 (in which case that
#' value is recycled) determining the start locations on the chromosome of the various guides.
#' @param num_nt_guides : an int determining the number of different NT guides to have. It can be 0 in which case
#' no NT guides are added.
#'
#' @return A \code{data.frame} with columns \code{grna_id}, \code{grna_target}, \code{chr}, \code{start}, and
#' \code{end} as described above. The values of  \code{grna_id}, \code{grna_target}, and \code{chr} indicate
#' how they relate. For instance, a \code{grna_id} of \code{"g2_t1_c3_d2"} is the 2nd guide for the target
#' \code{"t1_c3_d2"}; that target is the first target on the chromosome \code{"c3_d2"}, which used
#' \code{chr_distances[2]} and \code{chr_starts[2]} (unless the latter is a scalar in which case it's just
#' \code{chr_starts}).
#'
#' @examples
make_mock_grna_target_data <- function(num_guides_per_target, chr_distances, chr_starts, num_nt_guides) {
  if(length(chr_starts) == 1) {
    chr_starts <- rep(chr_starts, length(chr_distances))
  }
  if(length(chr_starts) != length(chr_distances)) {
    stop("`chr_starts` must be length 1 or have the same length as `chr_distances`.")
  }

  # first we make grna target data that is just doing different targets for different chromosomes
  # This data is formed by first having every target on a different chromosome,
  # and then repeating the guide counts but all on the same chromosome.
  target_and_chr_data <- .make_mock_target_and_chr_only(num_guides_per_target)

  # now we cross `target_and_chr_data` with each of the values in `chr_distance`
  # the result has the length of the final data and is just missing the grna_id column
  grna_target_data <- .make_mock_cross_target_chr_with_distances(
    target_and_chr_data = target_and_chr_data,
    chr_distances = chr_distances, chr_starts = chr_starts
  )

  # adding in unique names for each grna
  # this gets the unique target names, ensuring they stay in the order that they appear in the data
  target_labels <- with(grna_target_data, grna_target[!duplicated(grna_target)])
  # target_labels <- grna_target_data$grna_target[!duplicated(grna_target_data$grna_target)]
  # This vector enumerates each unique target: we get a sequence 1:(# times that target appears)
  grna_id_number <- table(grna_target_data$grna_target)[target_labels] |> sapply(function(i) 1:i) |> unlist()
  grna_target_data <- cbind(
    grna_id = paste0("g", grna_id_number, "_", grna_target_data$grna_target),
    grna_target_data
  )
  if(num_nt_guides > 0) {
    grna_target_data <- grna_target_data |>
      rbind(
        data.frame(
          grna_id = paste0("nt", 1:num_nt_guides),
          grna_target = "nt",
          chr = NA_character_, start = NA, end = NA
        )
      )
  }
  grna_target_data
}

######################################################################################################
##                                          Unit tests                                              ##
######################################################################################################

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

  for(input in inputs) {
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

  for(target_and_chr_input in target_and_chr_data_inputs) {
    for(chr_distance_and_start_input in chr_distance_and_start_inputs) {
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

  for(target_and_chr_input in target_and_chr_data_inputs) {
    for(chr_distance_and_start_input in chr_distance_and_start_inputs) {
      for(num_nt in num_nts) {
        out <- make_mock_grna_target_data(
          num_guides_per_target = target_and_chr_input,
          chr_distances = chr_distance_and_start_input$distance,
          chr_starts = chr_distance_and_start_input$start,
          num_nt_guides = num_nt
        )

        # `grna_target`, `chr`, `start` and `end` aren't modified
        # so this just tests `grna_id` and the appended NT guides
        out_pieces <- split(out, ifelse(out$grna_target == "nt", "nt", "non_nt"))

        if(num_nt > 0) {
          expect_equal(nrow(out_pieces$nt), num_nt)
          expect_equal(sub("nt", "", out_pieces$nt$grna_id) |> as.numeric(), 1:num_nt)
        }
        expect_equal(sub("g.*_(t)", "\\1", out_pieces$non_nt$grna_id),
                     out_pieces$non_nt$grna_target)
        with(out_pieces$non_nt, {
          for(target in unique(grna_target)) {
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

## example
# num_g_1 <- c(1,2,5)
# chr_dist_1 <- c(0, 1, 2, 100)
# num_nt_1 <- 2
#
# out1 <- make_mock_grna_target_data(
#   num_guides_per_target = num_g_1, chr_distances = chr_dist_1,
#   chr_starts = 1, num_nt_guides = num_nt_1
# )














