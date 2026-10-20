test_that("replacing guide counts refreshes summaries and invalidates prepared state atomically", {
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    case <- readRDS(test_path(
        "fixtures",
        "refactor",
        "cases.rds"
    ))$high_complement_crt
    object <- suppressMessages(refactor_prepare(inputs, case))
    object <- withr::with_seed(
        1901L,
        suppressMessages(run_discovery_analysis(
            object,
            print_progress = FALSE,
            parallel = FALSE
        ))
    )
    # Nonempty real fits/result make the invalidation checks non-vacuous.
    expect_gt(length(object@fits@response_fits), 0L)
    expect_gt(length(object@fits@perturbation_fits), 0L)
    expect_true(result_record_is_current(object, "run_discovery_analysis"))
    before <- serialize(object, NULL)
    original <- as.matrix(object@data@grna_matrix)
    replacement <- original
    replacement[1L, 1L] <- max(replacement[, 1L]) + 50
    changed <- set_grna_matrix(object, replacement)
    expect_identical(serialize(object, NULL), before)
    expect_equal(as.matrix(changed@data@grna_matrix), replacement)
    expect_equal(
        changed@data@covariate_data_frame$grna_n_umis,
        unname(colSums(replacement))
    )
    expect_equal(
        changed@data@covariate_data_frame$grna_n_nonzero,
        unname(colSums(replacement > .5))
    )
    expect_identical(
        changed@data@import_grna_assignment_info$max_grna[1L],
        rownames(replacement)[1L]
    )
    expect_equal(
        changed@data@import_grna_assignment_info$max_grna_frac_umis,
        unname(apply(replacement, 2L, max) / colSums(replacement))
    )
    expect_false(workflow_status(changed)[["assign_grnas"]])
    expect_false(workflow_status(changed)[["run_qc"]])
    expect_false(result_record_is_current(changed, "run_discovery_analysis"))
    expect_length(changed@fits@response_fits, 0L)
    expect_length(changed@fits@perturbation_fits, 0L)
    expect_identical(changed@data@response_matrix, object@data@response_matrix)
    response_replacement <- as.matrix(object@data@response_matrix)
    response_replacement[1L, 1L] <- response_replacement[1L, 1L] + 50
    changed_response <- set_response_matrix(object, response_replacement)
    expect_identical(serialize(object, NULL), before)
    expect_equal(
        as.matrix(changed_response@data@response_matrix),
        response_replacement
    )
    expect_equal(
        changed_response@data@covariate_data_frame$response_n_umis,
        unname(colSums(response_replacement))
    )
    expect_false(workflow_status(changed_response)[["run_qc"]])
    expect_false(result_record_is_current(
        changed_response,
        "run_discovery_analysis"
    ))
    expect_length(changed_response@fits@response_fits, 0L)
    expect_length(changed_response@fits@perturbation_fits, 0L)
    expect_identical(changed_response@data@grna_matrix, object@data@grna_matrix)
    for (modality in c("response", "grna")) {
        value <- as.matrix(methods::slot(
            object@data,
            paste0(modality, "_matrix")
        ))
        for (axis in c("features", "cells")) {
            bad <- value
            if (axis == "features") {
                rownames(bad)[1L] <- "different-feature"
            } else {
                colnames(bad)[1L] <- "different-cell"
            }
            setter <- if (modality == "grna") {
                set_grna_matrix
            } else {
                set_response_matrix
            }
            expect_error(
                setter(object, bad),
                "preserve aligned cells and feature identities"
            )
            expect_identical(serialize(object, NULL), before)
        }
    }
})

test_that("fit merges reject duplicate identities and explicitly remove invalidated groups", {
    cache <- list(A = list(coef = 1), B = list(coef = 2))
    before <- serialize(cache, NULL)
    expect_error(
        merge_response_fits(cache, list(A = list(coef = 3))),
        "must be unique"
    )
    expect_error(
        merge_response_fits(list(), list(list(coef = 1))),
        "must be unique"
    )
    for (updates in list(
        list(list(coef = 1)),
        setNames(list(1, 2), c("A", "A")),
        setNames(list(1), NA_character_),
        setNames(list(1), "")
    )) {
        expect_error(
            merge_perturbation_fits(cache, updates),
            "named and unique"
        )
    }
    actual <- merge_perturbation_fits(cache, list(A = NULL, C = list(coef = 3)))
    expect_identical(actual, list(B = list(coef = 2), C = list(coef = 3)))
    expect_identical(serialize(cache, NULL), before)
})
