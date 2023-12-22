set_response_matrix <- function(sceptre_object, response_matrix) {
  sceptre_object@response_matrix[[1]] <- response_matrix
  return(sceptre_object)
}

get_response_matrix <- function(sceptre_object) {
  return(sceptre_object@response_matrix[[1]] )
}

set_grna_matrix <- function(sceptre_object, grna_matrix) {
  sceptre_object@grna_matrix[[1]] <- grna_matrix
  return(sceptre_object)
}

get_grna_matrix <- function(sceptre_object) {
  return(sceptre_object@grna_matrix)
}
