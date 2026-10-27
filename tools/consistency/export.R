# Frozen behavior suite: language-neutral export (TSV and Matrix Market).
#
# Sourced together with harness.R and fixture-io.R. Writes only below the
# directory it is given. The format is described in README.md ("Fixture
# format"). The .rds observations remain the complete record; these files
# are a flattened copy for other languages, such as pysceptre's tests.

# Columns that hold resampled null statistics (output_amount = 3). They are
# recorded in the .rds observation only, to keep the TSV tables small.
consistency_tsv_omitted <- function(columns) grepl("^z_null_[0-9]+$", columns)

consistency_column_type <- function(x) {
    if (is.logical(x)) {
        "logical"
    } else if (is.integer(x)) {
        "integer"
    } else if (is.double(x)) {
        "double"
    } else if (is.character(x)) {
        "character"
    } else {
        stop("Unsupported column type: ", paste(class(x), collapse = "/"))
    }
}

consistency_format_column <- function(x) {
    out <- if (is.double(x)) {
        sprintf("%.17g", x)
    } else if (is.character(x)) {
        if (any(grepl("[\t\r\n]", x)) || any(x %in% "NA")) {
            stop("Character values may not contain tabs, newlines or 'NA'.")
        }
        x
    } else {
        as.character(x)
    }
    out[is.na(x)] <- "NA"
    out
}

# Write one table; return its manifest entry (types, rows, omitted columns).
consistency_write_tsv <- function(table, path) {
    omitted <- names(table)[consistency_tsv_omitted(names(table))]
    table <- table[, setdiff(names(table), omitted), drop = FALSE]
    types <- lapply(table, consistency_column_type)
    lines <- if (ncol(table) == 0L) {
        character()
    } else {
        do.call(
            paste,
            c(lapply(table, consistency_format_column), sep = "\t")
        )
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    connection <- file(path, open = "wb")
    on.exit(close(connection), add = TRUE)
    writeLines(
        enc2utf8(c(paste(names(table), collapse = "\t"), lines)),
        connection,
        sep = "\n",
        useBytes = TRUE
    )
    list(
        columns = types,
        n_rows = nrow(table),
        omitted_columns = if (length(omitted)) {
            list(pattern = "^z_null_[0-9]+$", count = length(omitted))
        }
    )
}

consistency_read_tsv <- function(path, columns) {
    lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
    header <- strsplit(lines[[1L]], "\t", fixed = TRUE)[[1L]]
    if (!identical(header, names(columns))) {
        stop("TSV header does not match the manifest: ", path)
    }
    fields <- strsplit(lines[-1L], "\t", fixed = TRUE)
    out <- lapply(seq_along(header), function(k) {
        values <- vapply(fields, `[[`, character(1), k)
        values[values == "NA"] <- NA_character_
        switch(
            columns[[k]],
            logical = as.logical(values),
            integer = as.integer(values),
            double = as.double(values),
            character = values
        )
    })
    names(out) <- header
    as.data.frame(out, stringsAsFactors = FALSE, optional = TRUE)
}

# A data frame for export: factors become character, list columns are dropped
# and non-automatic row names become a first column, `row_name`.
consistency_flat_table <- function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    row_names <- attr(x, "row.names")
    x <- x[, vapply(x, is.atomic, logical(1)), drop = FALSE]
    for (name in names(x)) {
        if (is.factor(x[[name]])) x[[name]] <- as.character(x[[name]])
    }
    if (is.character(row_names)) {
        x <- data.frame(row_name = row_names, x, check.names = FALSE)
    }
    rownames(x) <- NULL
    x
}

# Rebuild a data frame from consistency_plain() output, then flatten it.
consistency_unplain_table <- function(plain) {
    columns <- plain$columns
    out <- structure(
        columns,
        names = plain$names,
        row.names = plain$row_names,
        class = "data.frame"
    )
    consistency_flat_table(out)
}

# All `tables` lists inside an observation, keyed by their path.
consistency_find_tables <- function(value, path = character()) {
    if (!is.list(value) || is.data.frame(value)) {
        return(list())
    }
    found <- list()
    if (is.list(value$tables) && !is.data.frame(value$tables)) {
        for (name in names(value$tables)) {
            key <- paste(c(path, name), collapse = "__")
            found[[key]] <- value$tables[[name]]
        }
    }
    for (name in setdiff(names(value), "tables")) {
        found <- c(found, consistency_find_tables(value[[name]], c(path, name)))
    }
    found
}

# Write the case's tables; return their manifest entries.
consistency_export_case <- function(name, observation, out) {
    tables <- consistency_find_tables(observation$value)
    written <- list()
    for (key in names(tables)) {
        if (consistency_tsv_skipped(name, key)) {
            next
        }
        file <- file.path("tables", name, paste0(key, ".tsv"))
        written[[key]] <- c(
            list(file = file),
            consistency_write_tsv(
                consistency_unplain_table(tables[[key]]),
                file.path(out, file)
            )
        )
    }
    written
}

# Tables kept in the .rds observation only, to keep the TSV set small:
# count matrices of file imports (the inputs are standard files already),
# trans pairs (a full cross product) and the per-cell tables of nested
# observations (cache and edge cases repeat them for every analysis).
consistency_tsv_skipped <- function(case, key) {
    file_imports <- c(
        "import_cellranger_two_directories",
        "import_cellranger_one_directory",
        "import_parse"
    )
    (case %in% file_imports && key %in% c("response_matrix", "grna_matrix")) ||
        (case == "pairs_constructors" &&
            key %in% c("trans_none", "trans_pc_targets")) ||
        grepl("__(cell_covariates|grna_assignments)$", key)
}

consistency_write_dense_tsv <- function(m, path, id = "id") {
    table <- data.frame(rownames(m), m, check.names = FALSE)
    names(table)[[1L]] <- id
    consistency_write_tsv(table, path)
}

consistency_write_mtx_gz <- function(m, path) {
    m <- methods::as(methods::as(m, "CsparseMatrix"), "TsparseMatrix")
    connection <- gzfile(path, open = "wb")
    on.exit(close(connection), add = TRUE)
    is_integer <- all(m@x == round(m@x))
    writeLines(
        c(
            paste(
                "%%MatrixMarket matrix coordinate",
                if (is_integer) "integer" else "real",
                "general"
            ),
            paste(nrow(m), ncol(m), length(m@x)),
            paste(
                m@i + 1L,
                m@j + 1L,
                if (is_integer) {
                    format(m@x, scientific = FALSE, trim = TRUE)
                } else {
                    sprintf("%.17g", m@x)
                }
            )
        ),
        connection,
        sep = "\n"
    )
}

# Language-neutral copies of the inputs used by the cases.
consistency_write_inputs <- function(directory) {
    synthetic <- consistency_synthetic_inputs()
    s <- file.path(directory, "synthetic")
    consistency_write_dense_tsv(
        synthetic$response_matrix,
        file.path(s, "response_matrix.tsv"),
        "response_id"
    )
    for (moi in c("low", "high")) {
        consistency_write_dense_tsv(
            synthetic[[paste0(moi, "_grna_matrix")]],
            file.path(s, paste0(moi, "_grna_matrix.tsv")),
            "grna_id"
        )
    }
    covariates <- consistency_flat_table(synthetic$extra_covariates)
    names(covariates)[[1L]] <- "cell"
    consistency_write_tsv(covariates, file.path(s, "extra_covariates.tsv"))
    for (name in c(
        "grna_target_data_frame",
        "discovery_pairs",
        "positive_control_pairs"
    )) {
        consistency_write_tsv(
            consistency_flat_table(synthetic[[name]]),
            file.path(s, paste0(name, ".tsv"))
        )
    }
    consistency_write_tsv(
        data.frame(cell_index = synthetic$remove_cells),
        file.path(s, "remove_cells.tsv")
    )
    jsonlite::write_json(
        list(
            common = list(
                formula_object = "~ batch + z",
                assign_grnas = list(method = "thresholding", threshold = 40),
                run_qc = list(
                    response_n_umis_range = c(0, 1),
                    response_n_nonzero_range = c(0, 1),
                    n_nonzero_trt_thresh = 7L,
                    n_nonzero_cntrl_thresh = 7L,
                    additional_cells_to_remove = "remove_cells.tsv"
                ),
                run_calibration_check = list(
                    n_calibration_pairs = 4L,
                    calibration_group_size = 2L
                ),
                seed = consistency_case_seed,
                setting_arguments = c(
                    moi = "import_data(moi)",
                    control = "set_analysis_parameters(control_group)",
                    mechanism = "set_analysis_parameters(resampling_mechanism)",
                    integration = paste0(
                        "set_analysis_parameters(",
                        "grna_integration_strategy)"
                    ),
                    side = "set_analysis_parameters(side)",
                    approximation = paste0(
                        "set_analysis_parameters(",
                        "resampling_approximation)"
                    ),
                    output = "output_amount of each analysis"
                )
            ),
            settings = consistency_workflow_settings()
        ),
        file.path(s, "workflow_settings.json"),
        auto_unbox = TRUE,
        pretty = TRUE
    )

    example <- consistency_example_data()
    e <- file.path(directory, "example")
    dir.create(e, recursive = TRUE, showWarnings = FALSE)
    for (moi in c("low", "high")) {
        data <- example[[moi]]
        for (name in c("response_matrix", "grna_matrix")) {
            consistency_write_mtx_gz(
                data[[name]],
                file.path(e, paste0(moi, "moi_", name, ".mtx.gz"))
            )
            ids <- list(rownames(data[[name]]), colnames(data[[name]]))
            consistency_write_tsv(
                data.frame(id = ids[[1L]]),
                file.path(e, paste0(moi, "moi_", name, "_rows.tsv"))
            )
            consistency_write_tsv(
                data.frame(id = ids[[2L]]),
                file.path(e, paste0(moi, "moi_", name, "_cols.tsv"))
            )
        }
        consistency_write_tsv(
            consistency_flat_table(data$extra_covariates),
            file.path(e, paste0(moi, "moi_extra_covariates.tsv"))
        )
    }
    consistency_write_tsv(
        consistency_flat_table(example$low$grna_target_data_frame),
        file.path(e, "lowmoi_grna_target_data_frame.tsv")
    )
    consistency_write_tsv(
        consistency_flat_table(example$high_targets),
        file.path(e, "highmoi_grna_target_data_frame.tsv")
    )
    consistency_write_tsv(
        data.frame(gene_name = example$high$gene_names),
        file.path(e, "highmoi_gene_names.tsv")
    )
}
