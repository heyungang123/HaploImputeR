test_that("twinsGenerator validates seed_size", {
  hap_ref <- matrix(sample(0:1, 1000, replace = TRUE), nrow = 100)
  count_obj <- matrix(c(50, 50), nrow = 100, ncol = 2)
  
  expect_error(twinsGenerator(hap_ref, count_obj, seed_size = 200, batch_size = 2, prior_window_size = 50), 
               "cannot exceed")
  expect_error(twinsGenerator(hap_ref, count_obj, seed_size = -1, batch_size = 2, prior_window_size = 50), 
               "must be positive")
})

test_that("twinsGenerator validates batch_size", {
  hap_ref <- matrix(sample(0:1, 1000, replace = TRUE), nrow = 100)
  count_obj <- matrix(c(50, 50), nrow = 100, ncol = 2)
  
  expect_error(twinsGenerator(hap_ref, count_obj, seed_size = 10, batch_size = 0, prior_window_size = 50), 
               "must be positive")
  expect_error(twinsGenerator(hap_ref, count_obj, seed_size = 10, batch_size = 100, prior_window_size = 50), 
               "less than num_sites")
})

test_that("twinsGenerator reproducibility with seed", {
  hap_ref <- matrix(sample(0:1, 1000, replace = TRUE), nrow = 100)
  count_obj <- matrix(c(50, 50), nrow = 100, ncol = 2)
  
  result1 <- twinsGenerator(hap_ref, count_obj, seed_size = 10, batch_size = 2, 
                            prior_window_size = 50, seed = 42, verbose = FALSE)
  result2 <- twinsGenerator(hap_ref, count_obj, seed_size = 10, batch_size = 2, 
                            prior_window_size = 50, seed = 42, verbose = FALSE)
  expect_equal(result1, result2)
})

test_that("twinsGenerator returns correct dimensions", {
  hap_ref <- matrix(sample(0:1, 1000, replace = TRUE), nrow = 100)
  count_obj <- matrix(c(50, 50), nrow = 100, ncol = 2)
  
  result <- twinsGenerator(hap_ref, count_obj, seed_size = 10, batch_size = 2, 
                           prior_window_size = 50, verbose = FALSE)
  expect_equal(nrow(result), 100)
  expect_equal(ncol(result), 10)
})