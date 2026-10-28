test_that("construct_trans_pairs() leaves the caller's positive control table unchanged", {
    response <- matrix(
        1,
        3L,
        2L,
        dimnames = list(c("r1", "r2", "r3"), c("cell-a", "cell-b"))
    )
    metadata <- data.frame(
        grna_id = paste0("guide-", seq_len(4L)),
        grna_target = c("target-z", "pc-target", "target-y", "non-targeting")
    )
    object <- methods::new(
        "sceptre_object",
        response_matrix = list(response),
        grna_target_data_frame = metadata
    )
    targets <- c("target-z", "pc-target", "target-y")
    product <- do.call(
        rbind,
        lapply(rownames(response), function(response_id) {
            data.frame(grna_target = targets, response_id = response_id)
        })
    )
    pc_frame <- data.frame(
        grna_target = c("target-z", "target-y", "target-z"),
        response_id = c("r1", "r3", "r1")
    )
    object_before <- serialize(object, NULL)
    for (as_table in c(FALSE, TRUE)) {
        pc <- if (as_table) data.table::as.data.table(pc_frame) else pc_frame
        before <- serialize(pc, NULL)
        for (mode in c("none", "pc_pairs", "pairs_containing_pc_targets")) {
            expected <- product
            if (mode == "pc_pairs") {
                excluded <- paste(pc_frame$grna_target, pc_frame$response_id)
                expected <- expected[
                    !paste(expected$grna_target, expected$response_id) %in%
                        excluded,
                ]
            } else if (mode == "pairs_containing_pc_targets") {
                expected <- expected[
                    !expected$grna_target %in% pc_frame$grna_target,
                ]
            }
            rownames(expected) <- NULL
            actual <- construct_trans_pairs(object, pc, mode)
            expect_identical(class(actual), c("data.table", "data.frame"))
            expect_identical(lapply(actual, identity), as.list(expected))
            expect_identical(
                class(pc),
                class(if (as_table) data.table::data.table() else data.frame())
            )
            expect_identical(serialize(pc, NULL), before)
            expect_identical(serialize(object, NULL), object_before)
        }
    }
})
