test_that("assignment getter owns zero-based slots and preserves saved memberships", {
    for (assignments in list(list(single = c(2L, 5L)),
                             list(first = c(2L, 5L), second = c(1L, 4L)),
                             list(first = integer(), second = c(1L, 4L)))) {
        object <- methods::new("sceptre_object",
            data = methods::new("ScreenData",
                grna_matrix = list(matrix(0, length(assignments), 5L))),
            preparation = methods::new("PreparationRecord",
                initial_grna_assignment_list = assignments,
                cells_in_use = c(2L, 4L, 5L)),
            functs_called = c(assign_grnas = TRUE, run_qc = TRUE))
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
