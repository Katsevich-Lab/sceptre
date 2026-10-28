# Files for an Analysis (development plan, slide 5), readable from R and Python:
#
#   screen.h5mu           response and gRNA counts with per-cell data (MuData)
#   grna_targets.parquet  the gRNA table
#   state.h5ad            per-cell stage values in obs, the design matrix in
#                         obsm, and every other value in uns$state: one JSON
#                         text from jsonlite::serializeJSON(), which keeps R
#                         types, names, order and missing values exactly
#   tables/NAME.parquet   tables of the specification and stages
#   results/NAME.parquet  calibration, power and discovery results
#   fits.parquet          cached response model fits, one row per term
#
# Writing works on the slots of the equivalent sceptre_object, and reading
# rebuilds that object and converts it with as_analysis(), so the stage
# records and their hashes are made the same way in both directions.

#' Write an Analysis to a directory, or read it back
#'
#' @param analysis an `Analysis`
#' @param directory the directory to write to, which must be empty or not
#'   exist yet, or to read from
#' @return `write_analysis()` returns `directory`, invisibly; `read_analysis()`
#'   returns an `Analysis`
#' @keywords internal
write_analysis <- function(analysis, directory) {
    check_file_packages()
    # files left from an earlier analysis would be read back as part of this one
    if (length(list.files(directory, all.files = TRUE, no.. = TRUE)) > 0L) {
        stop("`directory` must be empty or not exist yet.", call. = FALSE)
    }
    for (folder in c("tables", "results")) {
        dir.create(
            file.path(directory, folder),
            recursive = TRUE,
            showWarnings = FALSE
        )
    }
    object <- as_sceptre_object(analysis)
    write_screen(analysis@screen, directory)
    n_cells <- ncol(analysis@screen)
    obs <- data.frame(row.names = cell_names(analysis@screen))
    obsm <- list()
    state <- list()
    for (name in state_slots()) {
        value <- methods::slot(object, name)
        if (is.data.frame(value)) {
            write_table(
                value,
                file.path(directory, "tables", paste0(name, ".parquet"))
            )
        } else if (name %in% cell_slots && is_cell_column(name, value)) {
            obs[[name]] <- if (name %in% cell_index_slots) {
                seq_len(n_cells) %in% value
            } else {
                value
            }
        } else if (
            identical(name, "covariate_matrix") && nrow(value) == n_cells
        ) {
            obsm$covariate_matrix <- unname(value)
            state$covariate_matrix_attributes <- list(
                terms = colnames(value),
                assign = attr(value, "assign"),
                contrasts = attr(value, "contrasts")
            )
        } else {
            state[[name]] <- formulas_as_text(value)
        }
    }
    strict(anndataR::write_h5ad(
        anndataR::AnnData(
            obs = obs,
            obsm = obsm,
            uns = list(state = jsonlite::serializeJSON(state, digits = I(17)))
        ),
        file.path(directory, "state.h5ad"),
        mode = "w"
    ))
    for (name in names(analysis@results)) {
        write_table(
            analysis@results[[name]],
            file.path(directory, "results", paste0(name, ".parquet"))
        )
    }
    if (length(analysis@fits) > 0L) {
        arrow::write_parquet(
            fits_table(analysis@fits),
            file.path(directory, "fits.parquet")
        )
    }
    write_checksums(directory)
    invisible(directory)
}

#' @rdname write_analysis
#' @keywords internal
read_analysis <- function(directory) {
    check_file_packages()
    check_checksums(directory)
    screen <- read_screen(directory)
    file <- strict(anndataR::read_h5ad(
        file.path(directory, "state.h5ad"),
        as = "InMemoryAnnData"
    ))
    state <- formulas_from_text(jsonlite::unserializeJSON(file$uns$state))
    object <- as_sceptre_object(do.call(
        methods::new,
        c(
            list("Analysis", screen = screen),
            analysis_parts(methods::new("sceptre_object"))
        )
    ))
    for (name in state_slots()) {
        initial <- methods::slot(object, name)
        value <- if (is.data.frame(initial)) {
            read_table(
                file.path(directory, "tables", paste0(name, ".parquet")),
                default = initial
            )
        } else if (name %in% names(file$obs)) {
            column <- file$obs[[name]]
            if (name %in% cell_index_slots) which(column) else column
        } else if (
            identical(name, "covariate_matrix") &&
                !is.null(file$obsm$covariate_matrix)
        ) {
            covariate_matrix(file, state$covariate_matrix_attributes)
        } else {
            state[[name]]
        }
        methods::slot(object, name, check = FALSE) <- value
    }
    for (name in names(analysis_layout$results)) {
        slot <- analysis_layout$results[[name]]
        path <- file.path(directory, "results", paste0(name, ".parquet"))
        methods::slot(object, slot, check = FALSE) <- read_table(
            path,
            default = methods::slot(object, slot)
        )
    }
    path <- file.path(directory, "fits.parquet")
    if (file.exists(path)) {
        object@response_precomputations <- fits_list(read_table(path))
    }
    as_analysis(object)
}

check_file_packages <- function() {
    rlang::check_installed(
        c("anndataR", "arrow", "jsonlite", "MuData", "MultiAssayExperiment"),
        "to write or read analysis files."
    )
}

# anndataR reports values it cannot write or read as warnings; treat them as
# errors, so that nothing is dropped silently.
strict <- function(code) {
    withCallingHandlers(code, warning = function(w) {
        stop(conditionMessage(w), call. = FALSE)
    })
}

# Slots stored in state.h5ad and tables/: all but the Screen, results and fits.
state_slots <- function() {
    setdiff(
        methods::slotNames("sceptre_object"),
        c(analysis_layout$screen, analysis_layout$results, analysis_layout$fits)
    )
}

# An MD5 checksum of every file, so that reading stops at a changed, missing
# or added file.
write_checksums <- function(directory) {
    files <- setdiff(
        list.files(directory, recursive = TRUE),
        "checksums.tsv"
    )
    sums <- tools::md5sum(file.path(directory, files))
    utils::write.table(
        data.frame(file = files, md5 = unname(sums)),
        file.path(directory, "checksums.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}

check_checksums <- function(directory) {
    expected <- utils::read.delim(
        file.path(directory, "checksums.tsv"),
        colClasses = "character"
    )
    actual <- unname(tools::md5sum(file.path(directory, expected$file)))
    present <- setdiff(list.files(directory, recursive = TRUE), "checksums.tsv")
    # a file added later under an expected name would otherwise be read
    changed <- c(
        expected$file[is.na(actual) | actual != expected$md5],
        setdiff(present, expected$file)
    )
    if (length(changed) > 0L) {
        stop(
            "Analysis files are missing, changed or added since they were ",
            "written: ",
            paste(changed, collapse = ", "),
            call. = FALSE
        )
    }
}

# Per-cell slots, stored as obs columns; index slots as logical columns, which
# hold sorted, distinct cell indices exactly. Anything else goes to the JSON.
cell_slots <- c(
    "grnas_per_cell",
    "cells_in_use",
    "cells_w_zero_or_twoplus_grnas"
)
cell_index_slots <- c("cells_in_use", "cells_w_zero_or_twoplus_grnas")

is_cell_column <- function(name, value) {
    length(value) > 0L &&
        (!name %in% cell_index_slots || !is.unsorted(value, strictly = TRUE))
}

cell_names <- function(screen) {
    names <- colnames(screen)
    if (is.null(names)) as.character(seq_len(ncol(screen))) else names
}

covariate_matrix <- function(file, attributes) {
    matrix <- file$obsm$covariate_matrix
    dimnames(matrix) <- list(rownames(file$obs), attributes$terms)
    attr(matrix, "assign") <- attributes$assign
    if (!is.null(attributes$contrasts)) {
        attr(matrix, "contrasts") <- attributes$contrasts
    }
    matrix
}

# Parquet keeps a data.frame's columns, types, missing values and row names,
# and a data.table's class; arrow returns other data.frames (and tibbles, which
# it does not tell apart) as tibbles. Tables without columns are still at their
# initial value and are not written.
write_table <- function(table, path) {
    if (ncol(table) > 0L) {
        arrow::write_parquet(table, path)
    }
}

read_table <- function(path, default = data.frame()) {
    if (!file.exists(path)) {
        return(default)
    }
    # mmap = FALSE: memory-mapped tables keep the file open, and Windows then
    # cannot overwrite or delete it while the analysis is in memory
    table <- arrow::read_parquet(path, mmap = FALSE)
    if (!inherits(table, "data.table")) {
        table <- as.data.frame(table)
    }
    table
}

# JSON cannot hold an R formula, so formulas are stored as their text.
formulas_as_text <- function(x) {
    if (inherits(x, "formula")) {
        text <- if (length(x) == 0L) {
            ""
        } else {
            paste(deparse(x, width.cutoff = 500L), collapse = " ")
        }
        return(structure(text, class = "formula_text"))
    }
    if (is.list(x) && !is.data.frame(x)) {
        x[] <- lapply(x, formulas_as_text)
    }
    x
}

formulas_from_text <- function(x) {
    if (inherits(x, "formula_text")) {
        text <- unclass(x)
        return(
            if (nzchar(text)) {
                stats::as.formula(text, env = globalenv())
            } else {
                methods::new("formula")
            }
        )
    }
    if (is.list(x) && !is.data.frame(x)) {
        x[] <- lapply(x, formulas_from_text)
    }
    x
}

# --- Screen -----------------------------------------------------------------

# Missing cell covariates: MuData 1.16 reads NA in a text column back as an
# empty string (a categorical column keeps it), so text columns are written as
# categoricals until that is fixed upstream; and h5mu stores a missing number as
# NaN, the AnnData convention, which is read back as NA.
write_screen <- function(screen, directory) {
    write_table(
        S4Vectors::metadata(screen)$grna_target_data_frame,
        file.path(directory, "grna_targets.parquet")
    )
    covariates <- as.data.frame(
        SummarizedExperiment::colData(screen),
        optional = TRUE
    )
    text <- names(covariates)[vapply(covariates, is.character, logical(1))]
    covariates[text] <- lapply(covariates[text], factor)
    SummarizedExperiment::colData(screen) <- S4Vectors::DataFrame(
        covariates,
        check.names = FALSE
    )
    meta <- S4Vectors::metadata(screen)
    meta$grna_target_data_frame <- NULL
    meta$covariate_text_columns <- text
    grna <- SingleCellExperiment::altExp(screen, "gRNA")
    SingleCellExperiment::altExps(screen) <- list()
    S4Vectors::metadata(screen) <- list()
    MuData::writeH5MU(
        MultiAssayExperiment::MultiAssayExperiment(
            list(response = screen, grna = grna),
            metadata = meta
        ),
        file.path(directory, "screen.h5mu")
    )
}

read_screen <- function(directory) {
    # MultiAssayExperiment warns that it rebuilds the sample map read from the
    # file; that is expected
    mae <- withCallingHandlers(
        MuData::readH5MU(file.path(directory, "screen.h5mu")),
        warning = function(w) {
            if (grepl("sampleMap", conditionMessage(w), fixed = TRUE)) {
                invokeRestart("muffleWarning")
            }
        }
    )
    # h5mu stores the counts as each modality's main matrix, read back unnamed;
    # sceptre keeps counts row-accessible (one response or gRNA per row)
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
    meta <- S4Vectors::metadata(mae)
    covariates <- as.data.frame(
        SummarizedExperiment::colData(screen),
        optional = TRUE
    )
    text <- as.vector(meta$covariate_text_columns)
    covariates[text] <- lapply(covariates[text], as.character)
    numeric <- vapply(covariates, is.double, logical(1))
    covariates[numeric] <- lapply(covariates[numeric], function(x) {
        replace(x, is.nan(x), NA)
    })
    SummarizedExperiment::colData(screen) <- S4Vectors::DataFrame(
        covariates,
        check.names = FALSE
    )
    S4Vectors::metadata(screen) <- list(
        low_moi = as.vector(meta$low_moi),
        grna_target_data_frame = read_table(file.path(
            directory,
            "grna_targets.parquet"
        )),
        covariate_row_names = as.vector(meta$covariate_row_names)
    )
    screen
}

# --- Fits -------------------------------------------------------------------

# One row per response and term; the negative binomial size (theta) is the row
# whose term is missing, so it cannot clash with a coefficient's name.
fits_table <- function(fits) {
    terms <- lapply(fits, function(fit) c(names(fit$fitted_coefs), NA))
    data.frame(
        response_id = factor(
            rep(names(fits), lengths(terms)),
            levels = names(fits)
        ),
        term = unlist(terms, use.names = FALSE),
        value = unlist(
            lapply(fits, function(fit) c(unname(fit$fitted_coefs), fit$theta)),
            use.names = FALSE
        )
    )
}

fits_list <- function(table) {
    lapply(split(table, table$response_id), function(rows) {
        is_theta <- is.na(rows$term)
        list(
            fitted_coefs = stats::setNames(
                rows$value[!is_theta],
                rows$term[!is_theta]
            ),
            theta = rows$value[is_theta]
        )
    })
}
