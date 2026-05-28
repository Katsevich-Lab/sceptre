# Deprecated function aliases kept for backward compatibility.

#' @keywords internal
#' @noRd
fisher_yates_samlper <- function(n_tot, M, B) {
  .Deprecated("fisher_yates_sampler")
  fisher_yates_sampler(n_tot = n_tot, M = M, B = B)
}
