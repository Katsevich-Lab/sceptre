# Files for an Analysis (development plan, slide 5), readable from R and Python:
#
#   screen.h5mu         the Screen: response and gRNA counts, cell covariates
#   state.h5ad          the specification and stage records: per-cell values in
#                       obs, the design matrix in obsm, the rest in uns
#   results/NAME.parquet  calibration, power and discovery result tables
#   fits.parquet        cached response model fits, one row per coefficient
#
# Stage records are stored in the forms AnnData defines (columns, tables,
# categoricals), not as serialized R objects, so pysceptre can read them.
# Slots that still hold their initial value are not stored.

#' Write an Analysis to a directory, or read it back
#'
#' @param analysis an `Analysis`
#' @param directory the directory to write to or read from
#' @return `write_analysis()` returns `directory`, invisibly; `read_analysis()`
#'   returns an `Analysis`
#' @keywords internal
write_analysis <- function(analysis, directory) {
    dir.create(
        file.path(directory, "results"),
        recursive = TRUE,
        showWarnings = FALSE
    )
    write_screen(analysis@screen, file.path(directory, "screen.h5mu"))
    anndataR::write_h5ad(
        state_anndata(analysis),
        file.path(directory, "state.h5ad"),
        mode = "w"
    )
    for (name in names(analysis@results)) {
        table <- analysis@results[[name]]
        if (ncol(table) > 0L) {
            arrow::write_parquet(
                table,
                file.path(directory, "results", paste0(name, ".parquet"))
            )
        }
    }
    if (length(analysis@fits) > 0L) {
        arrow::write_parquet(
            fits_table(analysis@fits),
            file.path(directory, "fits.parquet")
        )
    }
    invisible(directory)
}

#' @rdname write_analysis
#' @keywords internal
read_analysis <- function(directory) {
    screen <- read_screen(file.path(directory, "screen.h5mu"))
    state <- anndataR::read_h5ad(
        file.path(directory, "state.h5ad"),
        as = "InMemoryAnnData"
    )
    analysis <- do.call(
        methods::new,
        c(
            list("Analysis", screen = screen),
            analysis_parts(methods::new("sceptre_object"))
        )
    )
    analysis <- state_analysis(state, analysis)
    for (name in names(analysis@results)) {
        path <- file.path(directory, "results", paste0(name, ".parquet"))
        if (file.exists(path)) {
            analysis@results[[name]] <- as.data.frame(arrow::read_parquet(path))
        }
    }
    path <- file.path(directory, "fits.parquet")
    if (file.exists(path)) {
        analysis@fits <- fits_list(as.data.frame(arrow::read_parquet(path)))
    }
    analysis
}

# --- Screen -----------------------------------------------------------------

# MuData 1.16 reads missing values back changed: NA in a text column becomes
# the string "NA" (a categorical column keeps it), and NA in a numeric column
# becomes NaN. Until that is fixed upstream, text columns are written as
# categoricals, and both are restored on reading. (gRNA positions and cell
# covariates are never NaN, so NaN can be read as NA.)
text_as_categorical <- function(table) {
    text <- vapply(table, is.character, logical(1))
    table[text] <- lapply(table[text], factor)
    attr(table, "text_columns") <- names(table)[text]
    table
}

restore_missing <- function(table, text_columns) {
    table[text_columns] <- lapply(table[text_columns], as.character)
    numeric <- vapply(table, is.double, logical(1))
    table[numeric] <- lapply(table[numeric], function(x) {
        x[is.nan(x)] <- NA
        x
    })
    table
}

write_screen <- function(screen, path) {
    grna <- SingleCellExperiment::altExp(screen, "gRNA")
    SingleCellExperiment::altExps(screen) <- list()
    meta <- S4Vectors::metadata(screen)
    # h5mu sorts table columns by name, so keep the gRNA table's column order
    meta$grna_target_columns <- names(meta$grna_target_data_frame)
    meta$grna_target_data_frame <- text_as_categorical(
        meta$grna_target_data_frame
    )
    meta$grna_target_text_columns <- attr(
        meta$grna_target_data_frame,
        "text_columns"
    )
    covariates <- text_as_categorical(as.data.frame(
        SummarizedExperiment::colData(screen),
        optional = TRUE
    ))
    meta$covariate_text_columns <- attr(covariates, "text_columns")
    SummarizedExperiment::colData(screen) <- S4Vectors::DataFrame(
        covariates,
        check.names = FALSE
    )
    S4Vectors::metadata(screen) <- list()
    MuData::writeH5MU(
        MultiAssayExperiment::MultiAssayExperiment(
            list(response = screen, grna = grna),
            metadata = meta
        ),
        path
    )
}

read_screen <- function(path) {
    mae <- MuData::readH5MU(path)
    # sceptre keeps counts row-accessible (one response or gRNA per row)
    # (h5mu stores the counts as its main matrix, which is read back unnamed)
    row_counts <- function(experiment) {
        counts <- methods::as(
            SummarizedExperiment::assay(experiment, 1L),
            "RsparseMatrix"
        )
        SummarizedExperiment::assays(experiment) <- list(counts = counts)
        experiment
    }
    screen <- row_counts(mae[["response"]])
    SingleCellExperiment::altExp(screen, "gRNA") <- row_counts(mae[["grna"]])
    meta <- as_r_value(S4Vectors::metadata(mae))
    meta$grna_target_data_frame <- restore_missing(
        as.data.frame(meta$grna_target_data_frame[meta$grna_target_columns]),
        meta$grna_target_text_columns
    )
    SummarizedExperiment::colData(screen) <- S4Vectors::DataFrame(
        restore_missing(
            as.data.frame(
                SummarizedExperiment::colData(screen),
                optional = TRUE
            ),
            meta$covariate_text_columns
        ),
        check.names = FALSE
    )
    meta[c(
        "grna_target_columns",
        "grna_target_text_columns",
        "covariate_text_columns"
    )] <- NULL
    S4Vectors::metadata(screen) <- meta
    screen
}

# --- State: specification and stage records ---------------------------------

# A named list of cell-index vectors as a table with one row per (group, cell).
# The group is a categorical whose categories keep the list's order, empty
# groups included.
index_table <- function(x) {
    data.frame(
        group = factor(rep(names(x), lengths(x)), levels = names(x)),
        cell = unlist(x, use.names = FALSE)
    )
}

index_list <- function(table) {
    split(as.integer(table$cell), table$group)
}

# Assignments are lists whose elements are index lists or index vectors. h5ad
# sorts list elements by name, so the original order is stored with them.
assignment_tables <- function(x) {
    c(
        lapply(x, function(element) {
            if (is.list(element)) index_table(element) else element
        }),
        list(element_order = names(x))
    )
}

assignment_lists <- function(x) {
    lapply(x[as.vector(x$element_order)], function(element) {
        if (is.data.frame(element)) index_list(element) else as.vector(element)
    })
}

# A named vector as a two-column table.
named_table <- function(x) data.frame(name = names(x), value = unname(x))

named_vector <- function(table) stats::setNames(table$value, table$name)

# Each stored value with its writer and reader. Values still at their initial
# value are skipped on writing, and keep it on reading.
state_fields <- list(
    spec = list(
        formula_object = list(
            write = function(x) {
                paste(deparse(x, width.cutoff = 500L), collapse = " ")
            },
            read = function(x) stats::as.formula(x, env = globalenv())
        ),
        cellwise_qc_thresholds = list(write = identity, read = function(x) {
            as_r_value(x)[c(
                "response_n_umis_range",
                "response_n_nonzero_range",
                "p_mito_threshold"
            )]
        })
    ),
    assignment = list(
        initial_grna_assignment_list = list(
            write = index_table,
            read = index_list
        ),
        grna_assignments_raw = list(
            write = assignment_tables,
            read = assignment_lists
        ),
        grna_assignments = list(
            write = assignment_tables,
            read = assignment_lists
        )
    ),
    qc = list(
        cell_removal_metrics = list(write = named_table, read = named_vector)
    ),
    status = list(
        functs_called = list(write = named_table, read = named_vector)
    ),
    # h5ad text arrays cannot hold a missing value
    legacy = list(
        elements_to_analyze = list(
            write = function(x) {
                list(value = ifelse(is.na(x), "", x), missing = is.na(x))
            },
            read = function(x) {
                ifelse(as.vector(x$missing), NA_character_, as.vector(x$value))
            }
        )
    )
)

# Per-cell vectors, stored as obs columns. Index vectors (cells_in_use) are
# stored as a logical column marking the cells.
cell_fields <- c(
    grnas_per_cell = "assignment",
    cells_in_use = "qc",
    cells_w_zero_or_twoplus_grnas = "assignment"
)
cell_index_fields <- c("cells_in_use", "cells_w_zero_or_twoplus_grnas")

state_anndata <- function(analysis) {
    defaults <- analysis_parts(methods::new("sceptre_object"))
    n_cells <- ncol(analysis@screen)
    cell_names <- colnames(analysis@screen)
    if (is.null(cell_names)) {
        cell_names <- as.character(seq_len(n_cells))
    }
    groups <- c(
        list(spec = analysis@spec@settings),
        lapply(analysis@stages, methods::slot, name = "contents"),
        list(status = analysis@status, legacy = analysis@legacy)
    )
    default_groups <- c(
        list(spec = defaults$spec@settings),
        lapply(defaults$stages, methods::slot, name = "contents"),
        list(status = defaults$status, legacy = defaults$legacy)
    )
    obs <- data.frame(row.names = cell_names)
    obsm <- list()
    uns <- list()
    for (group in names(groups)) {
        for (name in names(groups[[group]])) {
            value <- groups[[group]][[name]]
            if (identical(value, default_groups[[group]][[name]])) {
                next
            }
            if (name %in% names(cell_fields)) {
                obs[[name]] <- if (name %in% cell_index_fields) {
                    seq_len(n_cells) %in% value
                } else {
                    value
                }
            } else if (identical(name, "covariate_matrix")) {
                obsm$covariate_matrix <- unname(value)
                uns$design$covariate_matrix_attributes <- list(
                    terms = colnames(value),
                    assign = attr(value, "assign"),
                    contrasts = attr(value, "contrasts")
                )
            } else {
                field <- state_fields[[group]][[name]]
                uns[[group]][[name]] <- if (is.null(field)) {
                    value
                } else {
                    field$write(value)
                }
            }
        }
    }
    uns$stage_hashes <- named_table(vapply(
        analysis@stages,
        methods::slot,
        character(1),
        name = "hash"
    ))
    anndataR::AnnData(obs = obs, obsm = obsm, uns = uns)
}

state_analysis <- function(state, analysis) {
    obs <- state$obs
    uns <- state$uns
    n_cells <- nrow(obs)
    set_value <- function(group, name, value) {
        if (group == "spec") {
            analysis@spec@settings[[name]] <<- value
        } else if (group %in% names(analysis@stages)) {
            analysis@stages[[group]]@contents[[name]] <<- value
        } else {
            methods::slot(analysis, group)[[name]] <<- value
        }
    }
    for (name in intersect(names(cell_fields), names(obs))) {
        value <- obs[[name]]
        if (name %in% cell_index_fields) {
            value <- which(value)
        }
        set_value(cell_fields[[name]], name, value)
    }
    if (!is.null(state$obsm$covariate_matrix)) {
        attrs <- uns$design$covariate_matrix_attributes
        matrix <- state$obsm$covariate_matrix
        dimnames(matrix) <- list(rownames(obs), as.vector(attrs$terms))
        attr(matrix, "assign") <- as.vector(attrs$assign)
        if (!is.null(attrs$contrasts)) {
            attr(matrix, "contrasts") <- attrs$contrasts
        }
        set_value("design", "covariate_matrix", matrix)
    }
    for (group in setdiff(names(uns), "stage_hashes")) {
        for (name in setdiff(
            names(uns[[group]]),
            "covariate_matrix_attributes"
        )) {
            field <- state_fields[[group]][[name]]
            value <- uns[[group]][[name]]
            set_value(
                group,
                name,
                if (is.null(field)) as_r_value(value) else field$read(value)
            )
        }
    }
    hashes <- named_vector(uns$stage_hashes)
    for (stage in names(analysis@stages)) {
        analysis@stages[[stage]]@hash <- hashes[[stage]]
    }
    analysis
}

# h5ad returns one-dimensional arrays for vectors; sceptre uses plain vectors.
as_r_value <- function(x) {
    if (is.array(x) && length(dim(x)) == 1L) {
        return(as.vector(x))
    }
    if (is.list(x) && !is.data.frame(x)) {
        return(lapply(x, as_r_value))
    }
    x
}

# --- Fits -------------------------------------------------------------------

# One row per response and term; "theta" is the negative binomial size.
fits_table <- function(fits) {
    rows <- lapply(fits, function(fit) {
        c(fit$fitted_coefs, theta = fit$theta)
    })
    data.frame(
        response_id = factor(
            rep(names(fits), lengths(rows)),
            levels = names(fits)
        ),
        term = unlist(lapply(rows, names), use.names = FALSE),
        value = unlist(rows, use.names = FALSE)
    )
}

fits_list <- function(table) {
    lapply(split(table, table$response_id), function(rows) {
        is_theta <- rows$term == "theta"
        list(
            fitted_coefs = stats::setNames(
                rows$value[!is_theta],
                rows$term[!is_theta]
            ),
            theta = rows$value[is_theta]
        )
    })
}
