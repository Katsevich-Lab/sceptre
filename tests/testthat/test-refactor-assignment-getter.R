test_that("assignment getter owns zero-based slots and preserves saved memberships", {
    for (assignments in list(list(single = c(2L, 5L)),
                             list(first = c(2L, 5L), second = c(1L, 4L)),
                             list(first = integer(), second = c(1L, 4L)))) {
        counts <- matrix(0, length(assignments), 5L,
            dimnames = list(names(assignments), NULL))
        for (i in seq_along(assignments)) counts[i, assignments[[i]]] <- 1
        metadata <- data.frame(grna_id = names(assignments),
            grna_target = names(assignments))
        cells <- c(2L, 4L, 5L)
        raw <- list(grna_group_idxs = assignments, indiv_nt_grna_idxs = list())
        retained <- list(grna_group_idxs = lapply(assignments, function(x) {
            match(intersect(x, cells), cells)
        }), indiv_nt_grna_idxs = list())
        object <- methods::new("sceptre_object",
            data = methods::new("ScreenData", identity = "assignment-getter-fixture",
                grna_target_data_frame = metadata,
                response_matrix = set_matrix_accessibility(matrix(1, 1L, 5L)),
                grna_matrix = set_matrix_accessibility(counts)),
            preparation = methods::new("PreparationRecord",
                guide_groups = names(assignments), design_matrix = matrix(1, 5L, 1L),
                guide_assignments = assignments,
                source_groups = raw, analysis_groups = retained,
                analysis_cells = cells))
        for (stage in c("design", "assignment", "cell_qc", "pair_qc")) {
            object <- refactor_private("commit_preparation_stage")(object, stage)
        }
        before <- serialize(object, NULL)
        expected <- matrix(FALSE, length(assignments), 5L,
                           dimnames = list(names(assignments), NULL))
        for (row in seq_along(assignments)) expected[row, assignments[[row]]] <- TRUE
        full <- get_grna_assignments(object)
        expect_s4_class(full, "ngRMatrix")
        expect_identical(as.matrix(full), expected)
        expect_identical(full@j, unlist(unname(assignments)) - 1L)
        expect_identical(full@p, c(0L, cumsum(lengths(unname(assignments)))))
        expect_identical(get_grna_assignments(object), full)
        # This getter historically drops the one-guide result to a logical
        # vector after QC. Keep that behavior rather than changing its API.
        qc <- get_grna_assignments(object, apply_cellwise_qc = TRUE)
        actual <- if (is.matrix(qc) || methods::is(qc, "Matrix")) as.matrix(qc) else qc
        expect_identical(actual, expected[, c(2L, 4L, 5L)])
        expect_identical(get_grna_assignments(object, apply_cellwise_qc = TRUE), qc)
        expect_identical(serialize(object, NULL), before)
    }
})
