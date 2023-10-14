B <- 50000
eps <- 0.01

test_that("Hybrid Fisher Yates/IWOR sampler", {
  N <- sample(seq(5, 15), 1)
  M <- sample(seq(4, 10), 1)
  m <- sample(seq(2, M - 1), 1)

  samples_cpp <- hybrid_fisher_iwor_sampler(N = N, m = m, M = M, B = B)
  samples <- samples_cpp |> synth_idx_list_to_matrix() |> t()
  
  rm(samples_cpp); gc()
  
  # criterion 1: uniqueness
  expect_true(all(apply(samples, 1, function(r) length(unique(r)) == M)))
  
  # criterion 2: the first m columns have minumum 0 and maxium N + m - 1
  expect_equal(min(samples[,seq(1:m)]), 0)
  expect_equal(max(samples[,seq(1:m)]), N + m - 1)
  
  # criterion 3: column i (for i > m) is bounded above by N + i - 1 and below by 0
  for (i in seq(m + 1, M)) {
    expect_lte(max(samples[,i]), N + i - 1)
    expect_gte(min(samples[,i]), 0)
  }
  
  # criterion 4: each entry in columns 1:m has probability m/N
  for (j in seq(0, N + m - 1)) {
    h <- mean(apply(X = samples[,1:m], 1, FUN = function(r) j %in% r))
    expect_lt(abs(mean(h) - m/(m+N)), eps)
  }
  
  # criterion 5: each entry in columns j = (m+1):M has probability j/(j + N)
  for (j in seq(m + 1, M)) {
    curr_samples <- samples[,1:j]
    for (i in seq(0, N + j - 1)) {
      h <- apply(X = curr_samples, MARGIN = 1, FUN = function(r) i %in% r)
      expect_lt(abs(mean(h) - j/(j + N)), eps)
    }
  }
})


test_that("Fisher-Yates sampler check", {
  n_tot <- sample(seq(10, 15), 1)
  M <- sample(seq(4, 8), 1)
  
  samples_cpp_2 <- fisher_yates_samlper(n_tot = n_tot, M = M, B = B)
  samples <- samples_cpp_2 |> synth_idx_list_to_matrix() |> t()
  rm(samples_cpp_2); gc()
  
  # criterion 0: sample length
  expect_equal(ncol(samples), M)
  
  # criterion 1: uniqueness
  expect_true(all(apply(samples, 1, function(r) length(unique(r)) == M)))
  
  # criterion 2: minumum 0, maximum M - 1
  expect_equal(min(samples), 0)
  expect_equal(max(samples), n_tot - 1)
  
  # criterion 3: equal probability of each element being in the sample
  for (j in seq(0, n_tot - 1)) {
    h <- mean(apply(X = samples, 1, FUN = function(r) j %in% r))
    expect_lt(abs(mean(h) - M/n_tot), eps)
  }
})
