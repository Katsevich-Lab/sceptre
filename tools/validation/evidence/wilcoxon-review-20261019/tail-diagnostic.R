# Lightweight independent actual-count/rank arithmetic, no sceptre installation.
# Reconstruct W08 seed821 fixed-library effect example exactly.
set.seed(821L)
n <- 400L
trt <- rep(1:4, length.out = n) <= 2L
batch <- stats::rbinom(n, 1, .5)
counts <- vapply(seq_len(n), function(i) {
    p <- if (trt[i]) .04 else .02
    stats::rmultinom(1L, 1000L, c(.01, p, .02, 1 - .01 - p - .02))[, 1L]
}, integer(4))
value <- log1p(10 * counts[2L, ])
x <- value[trt]
y <- value[!trt]
ranks <- rank(c(x, y))
nx <- length(x)
ny <- length(y)
N <- nx + ny
ties <- as.numeric(table(c(x, y)))
U <- sum(ranks[seq_len(nx)]) - nx * (nx + 1) / 2
variance <- nx * ny / 12 * ((N + 1) - sum(ties^3 - ties) / (N * (N - 1)))
z <- (U - nx * ny / 2 - sign(U - nx * ny / 2) * .5) / sqrt(variance)
base <- stats::wilcox.test(x, y, exact = FALSE, correct = TRUE)$p.value
stable <- 2 * stats::pnorm(-abs(z))
cat('Runtime: ', R.version.string, '\nBase two-sided p: ', base,
    '\nTie-corrected continuity-adjusted z: ', z,
    '\nStable two-sided normal tail: ', stable, '\n', sep = '')
# A zero result is host/version dependent; stable normal arithmetic need not be.
stopifnot(is.finite(z), is.finite(stable), stable > 0, stable < 1e-50)
