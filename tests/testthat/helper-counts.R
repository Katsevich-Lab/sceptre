# The response counts of a sceptre_object as a numeric matrix. Disk-backed
# counts are stored, and read, as integers; in-memory counts are doubles.
response_values <- function(sceptre_object) {
    out <- as.matrix(get_response_matrix(sceptre_object))
    storage.mode(out) <- "double"
    out
}
