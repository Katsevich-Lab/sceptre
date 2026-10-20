test_that("small ondisc stores choose valid chunks without changing counts or IDs", {
    skip_if_not_installed("ondisc")
    for (entries in c(2L, 1000L, 1001L)) {
        directory <- tempfile("small-backed-")
        dir.create(directory)
        on.exit(unlink(directory, recursive = TRUE), add = TRUE)
        response <- methods::new(
            "dgRMatrix",
            Dim = c(1L, entries + 1L),
            Dimnames = list("r", NULL),
            p = c(0L, entries),
            j = seq_len(entries) - 1L,
            x = as.double(seq_len(entries))
        )
        guide <- response
        rownames(guide) <- "g"
        before <- digest::digest(list(response, guide), algo = "sha256")
        backed <- storage_create_backed_pair(response, guide, directory, 19L)
        expect_identical(backed$response_matrix@integer_id, 19L)
        expect_identical(backed$grna_matrix@integer_id, 19L)
        expect_identical(
            as.numeric(backed$response_matrix[1, ]),
            c(as.double(seq_len(entries)), 0)
        )
        expect_identical(
            as.numeric(backed$grna_matrix[1, ]),
            c(as.double(seq_len(entries)), 0)
        )
        expect_identical(dim(backed$response_matrix), dim(response))
        expect_identical(
            before,
            digest::digest(list(response, guide), algo = "sha256")
        )
    }
})

test_that("unsupported empty or one-entry stores fail before either file is written", {
    skip_if_not_installed("ondisc")
    for (entries in 0:1) {
        directory <- tempfile("small-backed-reject-")
        dir.create(directory)
        on.exit(unlink(directory, recursive = TRUE), add = TRUE)
        response <- methods::new(
            "dgRMatrix",
            Dim = c(1L, 4L),
            Dimnames = list("r", NULL),
            p = c(0L, 2L),
            j = 0:1,
            x = c(1, 2)
        )
        guide <- methods::new(
            "dgRMatrix",
            Dim = c(1L, 4L),
            Dimnames = list("g", NULL),
            p = c(0L, as.integer(entries)),
            j = seq_len(entries) - 1L,
            x = rep(1, entries)
        )
        expect_error(
            storage_create_backed_pair(response, guide, directory, 19L),
            "fewer than two"
        )
        expect_length(list.files(directory, all.files = TRUE, no.. = TRUE), 0)
    }
})

test_that("public small backed import preserves counts across save and reopen", {
    skip_if_not_installed("ondisc")
    directory <- tempfile("small-public-backed-")
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    response <- matrix(
        c(1, 0, 2, 3, 0, 4),
        nrow = 2L,
        dimnames = list(c("r1", "r2"), c("c1", "c2", "c3"))
    )
    guide <- matrix(
        c(3, 0, 0, 2, 1, 0),
        nrow = 2L,
        dimnames = list(c("g1", "nt"), colnames(response))
    )
    targets <- data.frame(
        grna_id = rownames(guide),
        grna_target = c("A", "non-targeting")
    )
    before <- serialize(list(response, guide, targets), NULL)
    object <- suppressMessages(import_data(
        response_matrix = response,
        grna_matrix = guide,
        grna_target_data_frame = targets,
        moi = "low",
        use_ondisc = TRUE,
        directory_to_write = directory
    ))
    expect_true(validObject(object))
    expect_identical(
        get_response_matrix(object)@integer_id,
        get_grna_matrix(object)@integer_id
    )
    expect_identical(serialize(list(response, guide, targets), NULL), before)
    write_ondisc_backed_sceptre_object(object, directory)
    restored <- read_ondisc_backed_sceptre_object(
        file.path(directory, "sceptre_object.rds"),
        file.path(directory, "response.odm"),
        file.path(directory, "grna.odm")
    )
    expect_true(validObject(restored))
    expect_identical(restored@data@identity, object@data@identity)
    for (index in seq_len(nrow(response))) {
        expect_identical(
            as.numeric(get_response_matrix(restored)[index, ]),
            as.numeric(response[index, ])
        )
        expect_identical(
            as.numeric(get_grna_matrix(restored)[index, ]),
            as.numeric(guide[index, ])
        )
    }
    expect_identical(
        rownames(get_response_matrix(restored)),
        rownames(response)
    )
    expect_identical(rownames(get_grna_matrix(restored)), rownames(guide))
})
