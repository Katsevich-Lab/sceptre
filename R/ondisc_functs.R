##########################
# READ AND WRITE FUNCTIONS
##########################
#' Write or read an `ondisc`-backed `sceptre_object` (deprecated)
#'
#' `write_ondisc_backed_sceptre_object()` and
#' `read_ondisc_backed_sceptre_object()` are deprecated; use
#' [write_sceptre_object()] and [read_sceptre_object()]. `sceptre` no longer
#' uses `ondisc`: a disk-backed `sceptre_object` keeps its response counts in
#' the MuData file `screen.h5mu` (see [import_data()]).
#' `write_ondisc_backed_sceptre_object()` calls [write_sceptre_object()], and
#' `read_ondisc_backed_sceptre_object()` reads the `sceptre_object` written
#' to the directory of `sceptre_object_fp` with [read_sceptre_object()], its
#' response counts left on disk. The `.odm` files of earlier versions of
#' `sceptre` cannot be read; re-import the data.
#'
#' @param sceptre_object a `sceptre_object`
#' @param directory_to_write the directory in which to write the
#' `sceptre_object`
#' @param sceptre_object_fp file path to the `sceptre_object.rds` file written
#' by `write_ondisc_backed_sceptre_object()` or [write_sceptre_object()]
#' @param response_odm_file_fp,grna_odm_file_fp ignored
#'
#' @return `write_ondisc_backed_sceptre_object()` returns `NULL` invisibly, and
#' `read_ondisc_backed_sceptre_object()` returns a disk-backed
#' `sceptre_object`
#' @export
#' @examples
#' data(lowmoi_example_data)
#' sceptre_object <- import_data(
#'   response_matrix = lowmoi_example_data$response_matrix,
#'   grna_matrix = lowmoi_example_data$grna_matrix,
#'   grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
#'   extra_covariates = lowmoi_example_data$extra_covariates,
#'   moi = "low"
#' )
#' directory <- file.path(tempdir(), "deprecated_example")
#' # instead of write_ondisc_backed_sceptre_object()
#' write_sceptre_object(sceptre_object, directory)
#' # instead of read_ondisc_backed_sceptre_object()
#' sceptre_object <- read_sceptre_object(directory, backed = TRUE)
read_ondisc_backed_sceptre_object <- function(
    sceptre_object_fp,
    response_odm_file_fp = NULL,
    grna_odm_file_fp = NULL
) {
    .Deprecated("read_sceptre_object")
    if (!is.null(response_odm_file_fp) || !is.null(grna_odm_file_fp)) {
        warning(
            "`response_odm_file_fp` and `grna_odm_file_fp` are ignored; ",
            "the counts are read from screen.h5mu.",
            call. = FALSE
        )
    }
    sceptre_object_fp <- prepare_file_path(
        sceptre_object_fp,
        "`sceptre_object_fp`"
    )
    directory <- dirname(sceptre_object_fp)
    if (!file.exists(file.path(directory, "screen.h5mu"))) {
        stop(
            "The directory ",
            directory,
            " contains no screen.h5mu file. sceptre no longer reads the ",
            ".odm files of `ondisc`-backed objects; re-import the data, for ",
            "example with `import_data(..., backed = TRUE)`."
        )
    }
    read_sceptre_object(directory, backed = TRUE)
}


#' @rdname read_ondisc_backed_sceptre_object
#' @export
write_ondisc_backed_sceptre_object <- function(
    sceptre_object,
    directory_to_write
) {
    .Deprecated("write_sceptre_object")
    write_sceptre_object(sceptre_object, directory_to_write)
}


###########################
# PIPELINE HELPER FUNCTIONS
###########################
get_id_vect <- function(v, pod_size) {
    n_elements <- length(v)
    breaks <- round(n_elements / pod_size)
    if (breaks >= 2) {
        as.integer(cut(seq(1, n_elements), breaks))
    } else {
        rep(1L, n_elements)
    }
}

write_vector <- function(vector, file_name) {
    file_con <- file(file_name)
    on.exit(close(file_con), add = TRUE)
    writeLines(as.character(vector), file_con)
}
