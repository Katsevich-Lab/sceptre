# Build small valid sceptre objects just for the public pair constructors.
# Specify expected pairs directly, without model fitting or guide assignment.
refactor_pair_object <- function() {
    response_ids <- c("left", "far", "right", "edge-left", "edge-right", "other", "large")
    response <- matrix(1, length(response_ids), 2L,
                       dimnames = list(response_ids, c("cell-a", "cell-b")))
    metadata <- data.frame(
        grna_id = paste0("guide-", seq_len(6L)),
        grna_target = c("target-z", "target-z", "non-targeting", "pc-target", "target-y", "target-wide"),
        chr = c(rep("chr1", 4L), "chr2", "chr3"),
        start = c(80L, 100L, 90L, 90L, 30L, 1800000000L),
        end = c(100L, 120L, 110L, 110L, 30L, 1800000020L)
    )
    methods::new("sceptre_object", data = methods::new("ScreenData",
                 response_matrix = set_matrix_accessibility(response), grna_target_data_frame = metadata),
        preparation = methods::new("PreparationRecord", guide_groups = metadata$grna_target))
}

test_that("cis native positions follow the inclusive predicate in supplied order", {
    select <- refactor_private("compute_genes_within_distance")
    expect_identical(select(100L, c(90L, 1000L, 110L), 20L), c(1L, 3L))
    expect_identical(select(100L, c(120L, 80L, 121L, 79L, 100L, 80L), 20L),
                     c(1L, 2L, 5L, 6L))
    expect_identical(select(100L, integer(), 20L), integer())
    expect_identical(select(100L, c(100L, 99L), -1L), integer())

    # Compute expected distances using wider numbers before subtracting. The
    # first and last positions differ by more than a C++ int can hold, even
    # though this is a tiny input; no huge matrix is needed to test overflow.
    positions <- c(.Machine$integer.max, -.Machine$integer.max, 0L, 1L,
                   1800000010L, -1800000010L, 80L, 120L)
    before <- serialize(positions, NULL)
    for (midpoint in c(0L, 100L, 1800000010L, -1800000010L, .Machine$integer.max)) {
        for (distance in c(-1L, 0L, 20L, 1800000000L, .Machine$integer.max)) {
            expected <- which(abs(as.double(positions) - as.double(midpoint)) <= as.double(distance))
            expect_identical(select(midpoint, positions, distance), expected)
        }
    }
    expect_identical(serialize(positions, NULL), before)
    expect_error(select(NA_integer_, 1L, 1L), "fit native integer coordinates", fixed = TRUE)
    expect_error(select(1L, 1L, NA_integer_), "fit native integer coordinates", fixed = TRUE)
    expect_error(select(1L, c(1L, NA_integer_), 1L), "cannot be missing", fixed = TRUE)
})

test_that("public cis construction preserves row order and widens midpoint addition", {
    withr::local_seed(271828L)
    object <- refactor_pair_object()
    positions <- data.frame(
        response_id = rownames(get_response_matrix(object)),
        chr = c(rep("chr1", 5L), "chr2", "chr3"),
        position = c(90L, 1000L, 110L, 80L, 120L, 30L, 1800000010L)
    )
    pc <- data.frame(grna_target = "pc-target", response_id = "left")
    before <- lapply(list(object, positions, pc), serialize, connection = NULL)
    seed <- .Random.seed
    actual <- construct_cis_pairs(object, pc, 20L, positions)
    expect_identical(actual, data.frame(
        grna_target = c(rep("target-z", 4L), "target-y", "target-wide"),
        response_id = c("left", "right", "edge-left", "edge-right", "other", "large")
    ))
    expect_identical(lapply(list(object, positions, pc), serialize, connection = NULL), before)
    expect_identical(.Random.seed, seed)
})

test_that("trans exclusion owns caller tables and retains product order and class", {
    object <- refactor_pair_object()
    targets <- c("target-z", "pc-target", "target-y", "target-wide")
    responses <- rownames(get_response_matrix(object))
    product <- do.call(rbind, lapply(responses, function(response) {
        data.frame(grna_target = targets, response_id = response)
    }))
    pc_frame <- data.frame(grna_target = c("target-z", "target-y", "target-z"),
                           response_id = c("left", "other", "left"))
    object_before <- serialize(object, NULL)
    for (as_table in c(FALSE, TRUE)) {
        pc <- if (as_table) data.table::as.data.table(pc_frame) else pc_frame
        before <- serialize(pc, NULL)
        for (mode in c("none", "pc_pairs", "pairs_containing_pc_targets")) {
            expected <- product
            if (mode == "pc_pairs") {
                excluded <- paste(pc_frame$grna_target, pc_frame$response_id)
                expected <- expected[!paste(expected$grna_target, expected$response_id) %in% excluded, ]
            } else if (mode == "pairs_containing_pc_targets") {
                expected <- expected[!expected$grna_target %in% pc_frame$grna_target, ]
            }
            rownames(expected) <- NULL
            actual <- construct_trans_pairs(object, pc, mode)
            expect_identical(class(actual), c("data.table", "data.frame"))
            expect_identical(lapply(actual, identity), as.list(expected))
            expect_identical(serialize(pc, NULL), before)
            expect_identical(serialize(object, NULL), object_before)
        }
    }
    expect_identical(lapply(construct_trans_pairs(object, pairs_to_exclude = "pc_pairs"), identity),
                     as.list(product))
})
