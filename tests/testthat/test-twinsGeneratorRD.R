test_that("twinsGeneratorRD validates prior_window_size vector", {
  hap_ref <- matrix(sample(0:1, 1000, replace = TRUE), nrow = 100)
  hap_seeds <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10)
  count_obj <- matrix(c(50, 50), nrow = 100, ncol = 2)
  
  # Wrong length
  win_wrong <- rep(200, 50)
  expect_error(twinsGeneratorRD(hap_ref, hap_seeds, count_obj, batch_size = 2, prior_window_size = win_wrong), 
               "length equal to num_sites")
  
  # Contains negative values
  win_neg <- c(rep(200, 50), rep(-1, 50))
  expect_error(twinsGeneratorRD(hap_ref, hap_seeds, count_obj, batch_size = 2, prior_window_size = win_neg), 
               "non-negative")
})

test_that("twinsGeneratorRD reproducibility with seed", {
  hap_ref <- matrix(sample(0:1, 500, replace = TRUE), nrow = 50)
  hap_seeds <- matrix(sample(0:1, 50, replace = TRUE), nrow = 5)
  count_obj <- matrix(c(5, 5), nrow = 50, ncol = 2)
  win_sizes <- winSizeR2(hap_ref, 0.05, 30, 10)
  
  result1 <- twinsGeneratorRD(hap_ref, hap_seeds, count_obj, batch_size = 2, 
                              prior_window_size = win_sizes, seed = 42, verbose = FALSE)
  result2 <- twinsGeneratorRD(hap_ref, hap_seeds, count_obj, batch_size = 2, 
                              prior_window_size = win_sizes, seed = 42, verbose = FALSE)
  expect_equal(result1, result2)
})

test_that("twinsGeneratorRD returns correct dimensions", {
  hap_ref <- matrix(sample(0:1, 500, replace = TRUE), nrow = 50)
  hap_seeds <- matrix(sample(0:1, 50, replace = TRUE), nrow = 5)
  count_obj <- matrix(c(5, 5), nrow = 50, ncol = 2)
  win_sizes <- winSizeR2(hap_ref, 0.05, 30, 10)
  
  result <- twinsGeneratorRD(hap_ref, hap_seeds, count_obj, batch_size = 2, 
                             prior_window_size = win_sizes, verbose = FALSE)
  expect_equal(nrow(result), 50)
})

test_that("twinsGeneratorRD handles zero window sizes", {
  hap_ref <- matrix(sample(0:1, 500, replace = TRUE), nrow = 50)
  hap_seeds <- matrix(sample(0:1, 50, replace = TRUE), nrow = 5)
  count_obj <- matrix(c(5, 5), nrow = 50, ncol = 2)
  
  # Window sizes with zeros (should be converted to 1)
  win_zeros <- c(0, rep(20, 49))
  result <- twinsGeneratorRD(hap_ref, hap_seeds, count_obj, batch_size = 2, 
                             prior_window_size = win_zeros, seed = 42, verbose = FALSE)
  expect_equal(nrow(result), 50)
})