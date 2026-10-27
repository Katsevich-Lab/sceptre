# The language-neutral TSV tables are faithful copies of the .rds tables.

test_that("TSV tables equal the flattened tables of the .rds observations", {
    manifest <- consistency_manifest()
    directory <- consistency_fixture_dir()
    for (name in names(manifest$cases)) {
        entry <- manifest$cases[[name]]
        observation <- readRDS(file.path(directory, entry$rds))
        tables <- consistency_find_tables(observation$value)
        for (key in names(entry$tables)) {
            spec <- entry$tables[[key]]
            columns <- unlist(spec$columns)
            if (is.null(columns)) {
                columns <- stats::setNames(character(), character())
            }
            actual <- consistency_read_tsv(
                file.path(directory, spec$file),
                columns
            )
            flat <- consistency_unplain_table(tables[[key]])
            omitted <- consistency_tsv_omitted(names(flat))
            expected <- flat[, !omitted, drop = FALSE]
            # %.17g identifies each double exactly; R's parser (unlike a
            # correctly rounded one, such as Python's float) can be 1 ulp off.
            expect_identical(
                consistency_differences(actual, expected, tolerance = 4e-16),
                character(),
                info = paste(name, key)
            )
            expect_identical(
                if (any(omitted)) as.integer(spec$omitted_columns$count),
                if (any(omitted)) sum(omitted),
                info = paste(name, key)
            )
        }
    }
})
