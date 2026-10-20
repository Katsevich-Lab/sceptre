# Compare the scheduling contract against the prior repeated-scan definition.
# Nonmonotone keys and interleaved failed QC distinguish row/order errors.
unit_index_reference <- function(context, plan, n_work_units, permutations) {
    field <- if (permutations) "response_idx" else "group_idx"
    keys <- unique(context$index$pair_rows[[field]])
    if (!length(keys)) {
        stop("Distributed execution requires eligible hypothesis rows.")
    }
    n <- min(as.integer(n_work_units), length(keys))
    chunks <- split(
        keys,
        as.integer(floor((seq_along(keys) - 1L) * n / length(keys))) + 1L
    )
    eligible <- which(plan@pairs$pass_qc)
    lapply(seq_along(chunks), function(i) {
        selected <- as.integer(chunks[[i]])
        rows <- unlist(
            lapply(selected, function(key) {
                eligible[which(context$index$pair_rows[[field]] == key)]
            }),
            use.names = FALSE
        )
        list(
            id = sprintf("unit-%06d", i),
            keys = selected,
            row_ids = as.integer(rows)
        )
    })
}

test_that("preindexed units retain exact scheduling and eligible row identity", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$low_nt_perm
    object <- refactor_capture(refactor_prepare(inputs, case))$value
    plan <- refactor_private("prepare_execution_plan")(
        object,
        "run_discovery_analysis",
        NULL,
        NULL
    )
    pair_keys <- c(51L, 2L, 51L, 18L, 2L, 4L, 18L, 51L, 4L)
    pass <- c(
        TRUE,
        FALSE,
        TRUE,
        TRUE,
        FALSE,
        TRUE,
        TRUE,
        TRUE,
        TRUE,
        FALSE,
        TRUE,
        TRUE
    )
    plan@pairs <- data.frame(pass_qc = pass)
    context <- list(
        index = list(
            pair_rows = data.frame(
                response_idx = pair_keys,
                group_idx = rev(pair_keys)
            )
        )
    )
    set.seed(1806)
    rng <- .Random.seed
    for (permutations in c(TRUE, FALSE)) {
        for (units in c(1L, 2L, 3L, 4L, 99L)) {
            actual <- refactor_private("execution_units")(
                context,
                plan,
                units,
                permutations
            )
            expect_identical(
                actual,
                unit_index_reference(context, plan, units, permutations)
            )
            ids <- unlist(lapply(actual, `[[`, "row_ids"), use.names = FALSE)
            expect_identical(sort(ids), which(pass))
            expect_false(anyDuplicated(ids) > 0L)
        }
    }
    expect_identical(.Random.seed, rng)
    context$index$pair_rows <- data.frame(
        response_idx = integer(),
        group_idx = integer()
    )
    plan@pairs <- data.frame(pass_qc = logical())
    expect_error(
        refactor_private("execution_units")(context, plan, 2L, TRUE),
        "eligible hypothesis"
    )
})
