# Small valid sceptre objects for the pair constructors; expected pairs are
# specified directly, without model fitting or gRNA assignment.
cis_pair_object <- function(response_ids, metadata) {
    response <- matrix(
        1,
        length(response_ids),
        2L,
        dimnames = list(response_ids, c("cell-a", "cell-b"))
    )
    methods::new(
        "sceptre_object",
        response_matrix = list(response),
        grna_target_data_frame = metadata
    )
}

test_that("cis gene selection checks every position in supplied order", {
    expect_identical(
        compute_genes_within_distance(100L, c(90L, 1000L, 110L), 20L),
        c(1L, 3L)
    )
    expect_identical(
        compute_genes_within_distance(
            100L,
            c(120L, 80L, 121L, 79L, 100L, 80L),
            20L
        ),
        c(1L, 2L, 5L, 6L)
    )
    expect_identical(
        compute_genes_within_distance(100L, integer(), 20L),
        integer()
    )
    expect_identical(
        compute_genes_within_distance(100L, c(100L, 99L), -1L),
        integer()
    )
})

test_that("construct_cis_pairs() finds genes listed after an out-of-range gene", {
    response_ids <- c(
        "left",
        "far",
        "right",
        "edge-left",
        "edge-right",
        "other"
    )
    metadata <- data.frame(
        grna_id = paste0("guide-", seq_len(5L)),
        grna_target = c(
            "target-z",
            "target-z",
            "non-targeting",
            "pc-target",
            "target-y"
        ),
        chr = c(rep("chr1", 4L), "chr2"),
        start = c(80L, 100L, 90L, 90L, 30L),
        end = c(100L, 120L, 110L, 110L, 30L)
    )
    object <- cis_pair_object(response_ids, metadata)
    positions <- data.frame(
        response_id = response_ids,
        chr = c(rep("chr1", 5L), "chr2"),
        position = c(90L, 1000L, 110L, 80L, 120L, 30L)
    )
    pc <- data.frame(grna_target = "pc-target", response_id = "left")
    actual <- construct_cis_pairs(object, pc, 20L, positions)
    expect_identical(
        actual,
        data.frame(
            grna_target = c(rep("target-z", 4L), "target-y"),
            response_id = c("left", "right", "edge-left", "edge-right", "other")
        )
    )
})

test_that("construct_cis_pairs() finds all nearby genes in the GRCh37 table", {
    # gene_position_data_frame_grch37 is not sorted by position within each
    # chromosome; previously, only one of these nearby genes was paired
    utils::data(
        "gene_position_data_frame_grch37",
        package = "sceptre",
        envir = environment()
    )
    positions <- gene_position_data_frame_grch37[
        gene_position_data_frame_grch37$chr == "chr1",
    ][seq_len(300L), ]
    metadata <- data.frame(
        grna_id = c("guide-1", "guide-2"),
        grna_target = c("enhancer", "non-targeting"),
        chr = "chr1",
        start = 739137L,
        end = 739137L
    )
    object <- cis_pair_object(positions$response_id, metadata)
    actual <- construct_cis_pairs(
        object,
        distance_threshold = 500000L,
        response_position_data_frame = positions
    )
    expected <- positions$response_id[
        abs(positions$position - 739137) <= 500000
    ]
    expect_gt(length(expected), 1L)
    expect_identical(actual$response_id, expected)
    expect_true(all(actual$grna_target == "enhancer"))
})

test_that("cis gene selection does not overflow for large coordinates", {
    # The first two positions differ by more than a C++ int can hold
    positions <- c(
        .Machine$integer.max,
        -.Machine$integer.max,
        0L,
        1L,
        1800000010L,
        -1800000010L,
        80L,
        120L
    )
    for (midpoint in c(
        0L,
        100L,
        1800000010L,
        -1800000010L,
        .Machine$integer.max
    )) {
        for (distance in c(-1L, 0L, 20L, 1800000000L, .Machine$integer.max)) {
            expected <- which(
                abs(as.double(positions) - as.double(midpoint)) <=
                    as.double(distance)
            )
            expect_identical(
                compute_genes_within_distance(midpoint, positions, distance),
                expected
            )
        }
    }
})

test_that("construct_cis_pairs() computes the midpoint of large coordinates", {
    response_ids <- c("near", "far")
    metadata <- data.frame(
        grna_id = c("guide-1", "guide-2"),
        grna_target = c("target-wide", "non-targeting"),
        chr = "chr3",
        start = 1800000000L,
        end = 1800000020L
    )
    object <- cis_pair_object(response_ids, metadata)
    positions <- data.frame(
        response_id = response_ids,
        chr = "chr3",
        position = c(1800000010L, 1L)
    )
    actual <- construct_cis_pairs(
        object,
        distance_threshold = 20L,
        response_position_data_frame = positions
    )
    expect_identical(
        actual,
        data.frame(grna_target = "target-wide", response_id = "near")
    )
})
