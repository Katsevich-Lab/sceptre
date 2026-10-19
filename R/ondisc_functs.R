##########################
# READ AND WRITE FUNCTIONS
##########################
#' Write or read an `ondisc`-backed `sceptre_object`
#'
#' `write_ondisc_backed_sceptre_object()` and
#' `read_ondisc_backed_sceptre_object()` enable the writing and reading of
#' `ondisc`-backed `sceptre_object`s, respectively. First,
#' `write_ondisc_backed_sceptre_object()` writes an `ondisc`-backed
#' `sceptre_object` to disk, creating a file `sceptre_object.rds` in the
#' specified directory. Next, `read_ondisc_backed_sceptre_object()` reads and
#' initializes a `sceptre_object` from a `sceptre_object.rds` file,
#' `response.odm` file, and `grna.odm` file stored on disk.
#'
#' @param sceptre_object a `sceptre_object`
#' @param directory_to_write the directory in which to write the
#' `sceptre_object.rds` file
#' @param sceptre_object_fp file path to a `sceptre_object.rds` file
#' @param response_odm_file_fp file path to the backing `.odm` file for the
#' response modality
#' @param grna_odm_file_fp file path to the backing `.odm` file for the gRNA
#' modality
#'
#' @return `write_ondisc_backed_sceptre_object()` returns NULL, and
#' `read_ondisc_backed_sceptre_object()` returns an `ondisc`-backed
#' `sceptre_object`
#' @export
#' @examples
#' data(lowmoi_example_data)
#' # 1. create ondisc-backed sceptre_object
#' sceptre_object <- import_data(
#'   response_matrix = lowmoi_example_data$response_matrix,
#'   grna_matrix = lowmoi_example_data$grna_matrix,
#'   grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
#'   extra_covariates = lowmoi_example_data$extra_covariates,
#'   moi = "low",
#'   use_ondisc = TRUE,
#'   directory_to_write = tempdir()
#' )
#'
#' # 2. write
#' write_ondisc_backed_sceptre_object(
#'   sceptre_object = sceptre_object,
#'   directory_to_write = tempdir()
#' )
#'
#' # 3. read
#' rm(sceptre_object)
#' sceptre_object <- read_ondisc_backed_sceptre_object(
#'   sceptre_object_fp = paste0(tempdir(), "/sceptre_object.rds"),
#'   response_odm_file_fp = paste0(tempdir(), "/response.odm"),
#'   grna_odm_file_fp = paste0(tempdir(), "/grna.odm")
#' )
read_ondisc_backed_sceptre_object <- function(
    sceptre_object_fp,
    response_odm_file_fp,
    grna_odm_file_fp
) {
    # read in objects
    sceptre_object_fp <- prepare_ondisc_file_path(
        sceptre_object_fp,
        "`sceptre_object_fp`"
    )
    response_odm_file_fp <- prepare_ondisc_file_path(
        response_odm_file_fp,
        "`response_odm_file_fp`"
    )
    grna_odm_file_fp <- prepare_ondisc_file_path(
        grna_odm_file_fp,
        "`grna_odm_file_fp`"
    )
    sceptre_object <- migrate_sceptre_object(readRDS(sceptre_object_fp))
    matrices <- storage_reopen_backed_pair(
        response_odm_file_fp,
        grna_odm_file_fp,
        sceptre_object@data@integer_id
    )
    if (length(sceptre_object@data@backend)) {
        validate_backed_descriptor(matrices$response_matrix, sceptre_object@data@backend$response)
        validate_backed_descriptor(matrices$grna_matrix, sceptre_object@data@backend$grna)
    } else if (is.null(sceptre_object@results@records$legacy)) {
        stop("Backing store has no anchored content descriptor; reimport the screen.")
    }
    # Reattachment is a process-local handle change after verified identity;
    # it must not be treated as a mutation of the scientific data.
    sceptre_object@data@response_matrix <- list(matrices$response_matrix)
    sceptre_object@data@grna_matrix <- list(matrices$grna_matrix)
    if (!length(sceptre_object@data@backend)) {
        sceptre_object <- initialize_screen_identity(sceptre_object)
    }
    sceptre_object@data@backend <- list(
        response = backed_store_descriptor(matrices$response_matrix),
        grna = backed_store_descriptor(matrices$grna_matrix))
    return(sceptre_object)
}


#' @rdname read_ondisc_backed_sceptre_object
#' @export
write_ondisc_backed_sceptre_object <- function(
    sceptre_object,
    directory_to_write
) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    validate_backing_for_save(sceptre_object)
    # Keep the imported descriptor identity. Only the location may change after
    # verified reattachment; never bless changed content by writing fresh hashes.
    # Persist descriptors, never process-local live HDF5 handles.
    sceptre_object@data@grna_matrix <- list()
    sceptre_object@data@response_matrix <- list()
    # check that dir exists
    directory_to_write <- prepare_directory_to_write(directory_to_write)
    saveRDS(sceptre_object, file.path(directory_to_write, "sceptre_object.rds"))
}


#' Label consecutive pipeline work items with approximate task numbers
#'
#' Nextflow supplies an ordered vector and a positive target number of items per
#' task (pod). Return an integer task label for each item, keeping item order.
#' The number of tasks is rounded, so task size is approximate. Empty input gives
#' integer(). No random splitting is used.
#' @noRd
get_id_vect <- function(v, pod_size) {
    n_elements <- length(v)
    breaks <- round(n_elements / pod_size)
    if (breaks >= 2) {
        as.integer(cut(seq(1, n_elements), breaks))
    } else {
        rep(1L, n_elements)
    }
}

#' Write pipeline items as unquoted text lines in their supplied order
#'
#' Convert values to text, keeping repeats, and replace the destination file.
#' Close the connection on error as well as success. Return invisible NULL.
#' Nextflow scripts call this helper directly.
#' @noRd
write_vector <- function(vector, file_name) {
    file_con <- file(file_name)
    on.exit(close(file_con), add = TRUE)
    writeLines(as.character(vector), file_con)
}
