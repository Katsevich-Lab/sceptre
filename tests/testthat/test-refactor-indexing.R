indexing_screen <- function() {
    counts <- matrix(0, 2L, 12L, dimnames = list(c("r", "other"), NULL))
    counts[1L, c(4L, 8L, 11L)] <- 1
    counts[2L, ] <- seq_len(12L)
    object <- methods::new("sceptre_object")
    store <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "RsparseMatrix")
    object <- refactor_private("set_response_matrix")(object, store)
    object@low_moi <- TRUE
    object@control_group_complement <- FALSE
    object@cells_in_use <- c(2L, 4L, 5L, 8L, 11L)
    object@grna_assignments_raw <- list(
        grna_group_idxs = list(A = c(5L, 4L)),
        indiv_nt_grna_idxs = list(ntA = 8L, ntB = c(2L, 11L))
    )
    list(object = object, store = store, counts = counts,
         design = cbind(intercept = 1, source_cell = seq_len(12L)))
}

test_that("noncontiguous cells preserve actual extraction and native QC coordinates", {
    screen <- indexing_screen()
    before <- serialize(screen, NULL)
    object <- refactor_private("update_grna_assignments_given_qc")(screen$object)
    groups <- object@grna_assignments
    discovery <- c(groups$grna_group_idxs$A, groups$all_nt_idxs)
    calibration <- groups$indiv_nt_grna_idxs$ntB
    expect_identical(groups$grna_group_idxs$A, c(2L, 3L))
    expect_identical(groups$all_nt_idxs, c(4L, 1L, 5L))
    expect_identical(calibration, c(2L, 3L))
    y <- refactor_private("load_row")(screen$store, "r")
    discovery_y <- y[object@cells_in_use][discovery]
    calibration_y <- y[object@cells_in_use][groups$all_nt_idxs]
    expect_identical(discovery_y, c(1, 0, 1, 0, 1))
    expect_identical(calibration_y, c(1, 0, 1))
    expect_identical(unname(screen$design[object@cells_in_use, , drop = FALSE][discovery, 2L]),
                     c(4, 5, 8, 2, 11))
    expect_identical(sum(discovery_y[1:2] != 0), 1L)
    expect_identical(sum(discovery_y[3:5] != 0), 2L)
    expect_identical(sum(calibration_y[calibration] != 0), 1L)
    expect_identical(sum(calibration_y[-calibration] != 0), 1L)
    native <- refactor_private("compute_nt_nonzero_matrix_and_n_ok_pairs_v3")(
        j = screen$store@j, p = screen$store@p, n_cells_orig = 12L,
        n_cells_sub = 5L, grna_group_idxs = groups$grna_group_idxs,
        indiv_nt_grna_idxs = groups$indiv_nt_grna_idxs,
        all_nt_idxs = groups$all_nt_idxs, to_analyze_response_idxs = 1L,
        to_analyze_grna_idxs = 1L, control_group_complement = FALSE,
        cells_in_use = object@cells_in_use
    )
    expect_identical(native$n_nonzero_trt, 1L)
    expect_identical(native$n_nonzero_cntrl, 2L)
    expect_identical(unname(native$n_nonzero_mat[, 1L]), c(1L, 1L))
    expect_identical(serialize(screen, NULL), before)
    expect_identical(refactor_private("update_grna_assignments_given_qc")(screen$object)@grna_assignments,
                     groups)
})

test_that("observed native positions are one based and synthetic positions are zero based", {
    bridge <- refactor_native_bridge(test_path("fixtures", "refactor", "native-bridge.cpp"))
    for (view in list(list(y = c(1, 0, 1, 0, 1), trt = c(1L, 2L)),
                      list(y = c(1, 0, 1), trt = c(2L, 3L)))) {
        n <- length(view$y)
        pieces <- refactor_private("compute_precomputation_pieces")(
            view$y, matrix(1, n, 1L), fitted_coefs = 0, theta = 2,
            full_test_stat = TRUE
        )
        # Identical observed/synthetic memberships exercise the two distinct
        # bases in the production entry point, without relying on a sampler.
        draws <- rep(list(view$trt - 1L), 3L)
        result <- refactor_private("run_low_level_test_full_v4")(
            y = view$y, mu = pieces$mu, a = pieces$a, w = pieces$w,
            D = pieces$D, trt_idxs = view$trt, n_trt = length(view$trt),
            use_all_cells = FALSE,
            synthetic_idxs = bridge$refactor_native_pointer(draws, n),
            B1 = 3L, B2 = 0L, B3 = 0L, fit_parametric_curve = FALSE,
            return_resampling_dist = TRUE, side_code = 1L
        )
        expected_score <- sum(pieces$a[view$trt]) /
            sqrt(sum(pieces$w[view$trt]) - sum(rowSums(pieces$D[, view$trt, drop = FALSE])^2))
        expect_equal(result$z_orig, expected_score, tolerance = 1e-12)
        expect_equal(result$resampling_dist, rep(expected_score, 3L), tolerance = 1e-12)
        expect_identical(result$p, 1)
    }
})

test_that("legacy calibration unions and overlapping NT remapping stay explicit", {
    get_indices <- refactor_private("get_idx_vector_calibration_check")
    groups <- list(a = c(4L, 1L), b = c(1L, 5L))
    expect_identical(get_indices("a&b", groups, TRUE)$trt_idxs, c(4L, 1L, 1L, 5L))
    expect_identical(get_indices("a&b", groups, FALSE)$trt_idxs, c(4L, 1L, 5L))
    # The persisted legacy schema uses cumulative group lengths even when the
    # NT union overlaps. Characterize it here rather than silently repairing it.
    mapped <- refactor_private("update_indiv_grna_assignments_for_nt_cells")(groups)
    expect_identical(mapped$all_nt_idxs, c(4L, 1L, 5L))
    expect_identical(mapped$indiv_nt_grna_idxs, list(a = 1:2, b = 3:4))
})

test_that("canonical contexts compose source, analysis and fit coordinates once", {
    screen <- indexing_screen()
    object <- refactor_private("update_grna_assignments_given_qc")(screen$object)
    pairs <- data.frame(response_id = c("other", "r", "other"),
                        grna_group = rep("A", 3L))
    before <- serialize(list(object, pairs), NULL)
    prepare <- refactor_private("prepare_index_context")
    view_for <- refactor_private("make_fit_view")
    index <- prepare(object@cells_in_use, 12L, object@grna_assignments, pairs,
        rownames(screen$store), FALSE, FALSE, TRUE)
    expect_identical(index$analysis_to_source, c(2L, 4L, 5L, 8L, 11L))
    expect_identical(index$nt_to_analysis, c(4L, 1L, 5L))
    expect_identical(index$nt_groups_analysis, list(ntA = 4L, ntB = c(1L, 5L)))
    view <- view_for(index, "A")
    expect_identical(view$source_cells, c(4L, 5L, 8L, 2L, 11L))
    expect_identical(view$trt_fit_pos, 1:2)
    expect_identical(view$control_fit_pos, 3:5)
    expect_identical(view$n_trt, 2L)
    expect_identical(unname(screen$counts[1L, view$source_cells]), c(1, 0, 1, 0, 1))
    expect_identical(screen$design[view$source_cells, , drop = FALSE],
        screen$design[c(4L, 5L, 8L, 2L, 11L), , drop = FALSE])

    # The request-row map survives traversal in response order, including
    # repeated row identities. Public duplicate rejection remains upstream.
    expect_identical(index$pair_rows$response_idx, c(2L, 1L, 2L))
    expect_identical(index$pair_rows$pair_row_id, 1:3)
    traversal <- unname(unlist(index$rows_by_response))
    expect_identical(traversal, c(2L, 1L, 3L))
    expect_identical(index$pair_rows$pair_row_id[traversal][order(traversal)], 1:3)

    calibration_pairs <- pairs
    calibration_pairs$grna_group <- "ntB"
    calibration <- prepare(object@cells_in_use, 12L, object@grna_assignments,
        calibration_pairs, rownames(screen$store), TRUE, FALSE, TRUE)
    nt_view <- view_for(calibration, "ntB")
    expect_identical(nt_view$source_cells, c(8L, 2L, 11L))
    expect_identical(nt_view$trt_fit_pos, c(2L, 3L))
    expect_identical(nt_view$source_cells[nt_view$control_fit_pos], 8L)
    expect_identical(refactor_private("prepared_index_factory")(calibration)("ntB"),
                     list(trt_idxs = c(2L, 3L), n_trt = 2L))
    expect_identical(view_for(index, "A"), view)
    expect_identical(serialize(list(object, pairs), NULL), before)
})

test_that("analysis-cell normalization preserves reordered requested cells", {
    make_map <- refactor_private("prepare_cell_map")
    normalize <- refactor_private("source_membership_to_analysis")
    map <- make_map(c(11L, 4L, 2L, 8L, 5L), 12L)
    expect_identical(map, c(11L, 4L, 2L, 8L, 5L))
    expect_identical(normalize(c(5L, 4L, 5L), map), c(2L, 5L))
    expect_identical(normalize(integer(), map), integer())
    expect_error(make_map(c(2L, 2L), 12L), "unique integer")
    expect_error(make_map(c(0L, 4L), 12L), "unique integer")
    expect_error(make_map(c(4, 4.5), 12L), "unique integer")
    expect_error(make_map(c(4L, NA_integer_), 12L), "unique integer")
    expect_identical(refactor_private("native_observed_positions")(c(2L, 5L)), c(2L, 5L))
})
