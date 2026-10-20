# Lightweight harness integrity checks; no package fit/native compilation.
source("tools/validation/scientific-harness.R")
cases <- validation_cases()
stopifnot(!anyDuplicated(cases$id))
x <- validation_input(cases[1L, ], 41L)
y <- validation_input(cases[1L, ], 41L)
z <- validation_input(cases[1L, ], 43L)
stopifnot(identical(x, y), !identical(x$response, z$response),
          identical(unname(x$truth), c(0, log(.5), log(2))),
          all(colSums(x$guides > 0) == 1L),
          all((x$guides[1L, ] > 0) == as.logical(x$assignment)))
high <- validation_input(cases[3L, ], 41L)
stopifnot(any(colSums(high$guides > 0) > 1L))
rows <- data.frame(case = "toy", response_id = "null", true_log_effect = 0,
                   error = c("", "", "", "", "failed"),
                   pass_qc = c(TRUE, TRUE, TRUE, FALSE, FALSE),
                   p_value = c(.01, .9, 2, NA, NA))
s <- validation_summary(rows, .05)
stopifnot(s$attempted == 5L, s$valid == 2L, s$invalid_p == 1L,
          s$errors == 1L, s$qc_failed == 1L, s$rejected == 1L,
          s$estimate == .5, s$rejection_per_attempt == .2,
          s$lower_95 < .5, s$upper_95 > .5)
stopifnot(all(is.na(validation_interval(0L, 0L))))
cat("Harness integrity checks passed.\n")
