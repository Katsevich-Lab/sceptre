run_reduced_em_algo_v2 <- function(pi_guesses, g_pert_guesses, g, g_mus_pert0) {
  # define options as well as "outer" Ti1, g_pert, and pi
  ep_tol <- 0.5 * 1e-4
  max_it <- 50
  min_it <- 3
  n <- length(g)

  outer_Ti1s <- numeric()
  outer_converged <- FALSE
  outer_g_pert <- NA_integer_
  outer_pi <- NA_integer_
  outer_log_lik <- -Inf
  outer_i <- NA_integer_

  # loop over initial guesses
  for (i in seq_along(pi_guesses)) {
    # define some useful variables
    converged <- FALSE
    prev_log_lik <- -Inf
    curr_g_pert <- g_pert_guesses[i]
    curr_pi <- pi_guesses[i]
    iteration <- 1L

    # Iterate through E and M steps until convergence.
    while (TRUE) {
      ##############
      # E step start
      ##############
      g_mus_pert1 <- g_mus_pert0 + curr_g_pert

      # compute log likelihood
      p0 <- exp(log(1 - curr_pi) + stats::dpois(g, g_mus_pert0, log = TRUE))
      p1 <- exp(log(curr_pi) + stats::dpois(g, g_mus_pert1, log = TRUE))
      s <- p0 + p1
      s[s < 1e-100] <- 1e-100
      curr_log_lik <- sum(log(s))

      # compute membership probabilities
      quotient <- log(1 - curr_pi) - log(curr_pi) + g * (log(g_mus_pert0) - log(g_mus_pert1)) + g_mus_pert1 - g_mus_pert0
      Ti1s <- 1/(exp(quotient) + 1)

      # verify that Ti1s are not NA or 0
      if (any(is.na(Ti1s)) || all(Ti1s == 0)) {
        curr_log_lik <- -Inf
        break()
      }

      # estimate new_pi and ensure new_pi is less than 0.5
      curr_pi <- sum(Ti1s)/n
      if (curr_pi > 0.5) {
        Ti1s <- 1 - Ti1s
        curr_pi <- 1 - curr_pi
      }

      # Check for convergence
      tol <- compute_tolerance(curr_log_lik, prev_log_lik)
      if (tol < ep_tol && iteration >= min_it) {
        converged <- TRUE
        break()
      } else {
        prev_log_lik <- curr_log_lik
        iteration <- iteration + 1L
        if (iteration >= max_it) break()
      }
      ############
      # E step end
      ############

      # M step
      e1 <- sum(Ti1s * g)
      e2 <- sum(Ti1s * g_mus_pert0)
      curr_g_pert <- log(e1) - log(e2)
    }

    if (curr_log_lik > outer_log_lik && converged) {
      outer_Ti1s <- Ti1s
      outer_i <- i
      outer_log_lik <- curr_log_lik
      outer_converged <- converged
    }
  }
  return(list(outer_i = outer_i, outer_log_lik = outer_log_lik,
              Ti1s = outer_Ti1s, outer_converged = outer_converged))
}


compute_tolerance <- function(curr_log_lik, prev_log_lik) {
  if (curr_log_lik == -Inf || prev_log_lik == -Inf) {
    tol <- Inf
  } else {
    tol <- abs(curr_log_lik - prev_log_lik)/min(abs(curr_log_lik), abs(prev_log_lik))
  }
  return(tol)
}
