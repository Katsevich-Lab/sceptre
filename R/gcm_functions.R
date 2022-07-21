#' Run GCM inference
#'
#' Run GCM inference using the gene and gRNA precomputations
#'
#' @param gene_resids the residuals of the gene regression
#' @param gRNA_resids the residuals of the gRNA regression
#' @param side sidedness of the test; one of "left", "right", or "both"
#'
#' @return a data frame containing columns "p_value" and "z_value"
run_sceptre_using_gcm <- function(gene_resids, gRNA_resids, side) {
  # obtain basic quantities
  r <- gene_resids * gRNA_resids
  n <- length(r)
  s_r <- sum(r)
  s_r2 <- sum(r^2)
  # compute the z-score (under the null hypothesis)
  top <- 1/sqrt(n) * s_r
  bottom <- sqrt(1/n * s_r2 - (1/n * s_r)^2)
  z <- top/bottom
  # compute the p-value
  p_val <- switch(side,
                  "left" = pnorm(q = z, lower.tail = TRUE),
                  "right" = pnorm(q = z, lower.tail = FALSE),
                  "both" = 2 * pnorm(q = abs(z), lower.tail = FALSE))
  return(data.frame(p_value = p_val, z_value = z))
}

