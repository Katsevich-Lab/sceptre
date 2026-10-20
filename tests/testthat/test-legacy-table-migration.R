legacy_table_fixture <- function(value) {
    # An unavailable S3 extension mimics a serialized tibble without depending
    # on any package's old-class registration in the test process.
    attr(value, "class") <- c(
        "sceptre_unregistered_table",
        "tbl_df",
        "tbl",
        "data.frame"
    )
    attr(value, "source_note") <- "preserve this attribute"
    value
}

legacy_metadata_schema_fixture <- function(version) {
    file <- if (version == 1L) {
        "sceptre-0.99.0-21f9ea0.rds"
    } else {
        "sceptre-schema2-4aad80b.rds"
    }
    readRDS(test_path("fixtures", "saved-objects", file))
}

test_that("legacy table normalization changes only the S3 container class", {
    original <- data.frame(
        grna_id = c("g2", "g1"),
        grna_target = c("B", "A"),
        integer_value = c(2L, NA_integer_),
        factor_value = factor(c("b", "a")),
        date_value = as.Date(c("2026-10-20", "2026-10-19")),
        check.names = FALSE,
        row.names = c("second", "first")
    )
    original$list_value <- list(list(id = "x"), 1:3)
    original <- legacy_table_fixture(original)
    before <- serialize(original, NULL)
    result <- legacy_slot_data_frame(original)
    expected <- original
    attr(expected, "class") <- "data.frame"
    expect_identical(result, expected)
    expect_identical(class(result), "data.frame")
    expect_identical(rownames(result), c("second", "first"))
    for (field in names(original)) {
        expect_identical(result[[field]], original[[field]])
    }
    expect_identical(serialize(original, NULL), before)
    wrapped <- asS4(original, complete = FALSE)
    wrapped_before <- serialize(wrapped, NULL)
    expect_true(isS4(wrapped))
    expect_identical(legacy_slot_data_frame(wrapped), expected)
    expect_false(isS4(legacy_slot_data_frame(wrapped)))
    expect_identical(serialize(wrapped, NULL), wrapped_before)
    expect_null(legacy_slot_data_frame(NULL))
    expect_identical(legacy_slot_data_frame(1:3), 1:3)
})

test_that("legacy table subclasses cross all typed component/request boundaries", {
    for (version in c(1L, 2L)) {
        plain <- legacy_metadata_schema_fixture(version)
        baseline <- migrate_sceptre_object(plain)
        old <- plain
        if (version == 1L) {
            for (field in c(
                "grna_target_data_frame",
                "covariate_data_frame",
                "discovery_pairs",
                "positive_control_pairs"
            )) {
                attr(old, field) <- asS4(
                    legacy_table_fixture(attr(old, field)),
                    complete = FALSE
                )
            }
        } else {
            for (owner in c("data", "spec", "preparation")) {
                component <- attr(old, owner)
                for (field in names(attributes(component))) {
                    value <- attr(component, field)
                    if (is.data.frame(value)) {
                        attr(component, field) <- asS4(
                            legacy_table_fixture(value),
                            complete = FALSE
                        )
                    }
                }
                if (owner == "spec") {
                    requested <- attr(component, "requested")
                    for (field in c(
                        "discovery_pairs",
                        "positive_control_pairs"
                    )) {
                        if (is.data.frame(requested[[field]])) {
                            requested[[field]] <- asS4(
                                legacy_table_fixture(requested[[field]]),
                                complete = FALSE
                            )
                        }
                    }
                    attr(component, "requested") <- requested
                }
                attr(old, owner) <- component
            }
        }
        before <- serialize(old, NULL)
        migrated <- migrate_sceptre_object(old)
        expect_true(validObject(migrated))
        for (field in c("grna_target_data_frame", "covariate_data_frame")) {
            table <- slot(migrated@data, field)
            reference <- slot(baseline@data, field)
            attr(table, "source_note") <- NULL
            expect_identical(table, reference)
        }
        for (owner in c("spec", "preparation")) {
            for (field in c("discovery_pairs", "positive_control_pairs")) {
                table <- slot(slot(migrated, owner), field)
                reference <- slot(slot(baseline, owner), field)
                attr(table, "source_note") <- NULL
                expect_identical(table, reference)
            }
        }
        expect_identical(serialize(old, NULL), before)
    }
})
