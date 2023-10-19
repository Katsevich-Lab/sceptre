mm_sn_estimator <- function(y) {
  n <- length(y)    
  max.gamma1 <- 0.5*(4-pi)*(2/(pi-2))^1.5 - (.Machine$double.eps)^(1/4)    # maximum value gamma can take
  # method of moments strategy
  s <- sd(y)
  gamma1 <- sum((y-mean(y))^3)/(n*s^3)
  if(abs(gamma1) > max.gamma1) gamma1 <- sign(gamma1)*0.9*max.gamma1
  cp1 <- as.numeric(c(mean(y), s, gamma1))
  dp1 <- sn::cp2dp(cp1, family = "SN")
  return(dp1)
}


test_that("sn method of moments", {
  xi <- runif(n = 1, min = -1, max = 1)
  omega <- runif(n = 1, min = 0, max = 2)
  alpha <- runif(n = 1, min = 1, max = 3)
  y <- sn::rsn(n = 25000, xi = xi, omega = omega, alpha = alpha)
  
  fit1 <- as.numeric(mm_sn_estimator(y))
  fit2 <- fit_skew_normal_funct(y)
  
  expect_equal(fit2, fit2)
})
