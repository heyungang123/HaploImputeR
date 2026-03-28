test_that("twinsGeneratorM handles different chromosome counts", {
  # This design allows using a larger reference population for model training
  # while generating a smaller target population
  hap_ref <- matrix(sample(0:1, 1000, replace = TRUE), nrow = 100, ncol = 10)
  hap_seeds <- matrix(sample(0:1, 50, replace = TRUE), nrow = 5, ncol = 10)
  count_obj <- matrix(c(5, 5), nrow = 100, ncol = 2)
  
  # Different column counts should work
  hap_ref_large <- matrix(sample(0:1, 2000, replace = TRUE), nrow = 100, ncol = 20)
  result <- twinsGeneratorM(hap_ref_large, hap_seeds, count_obj, batch_size = 2, 
                            prior_window_size = 50, seed = 42, verbose = FALSE)
  expect_equal(dim(result), c(100, 10))  # Result has same columns as hap_seeds
})

test_that("twinsGeneratorM reproducibility with seed", {
  hap_ref <- matrix(sample(0:1, 1000, replace = TRUE), nrow = 100)
  hap_seeds <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10)
  count_obj <- matrix(c(50, 50), nrow = 100, ncol = 2)
  
  result1 <- twinsGeneratorM(hap_ref, hap_seeds, count_obj, batch_size = 2, 
                             prior_window_size = 50, seed = 42, verbose = FALSE)
  result2 <- twinsGeneratorM(hap_ref, hap_seeds, count_obj, batch_size = 2, 
                             prior_window_size = 50, seed = 42, verbose = FALSE)
  expect_equal(result1, result2)
})

test_that("twinsGeneratorM returns correct dimensions", {
  hap_ref <- matrix(sample(0:1, 1000, replace = TRUE), nrow = 100)
  hap_seeds <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10)
  count_obj <- matrix(c(50, 50), nrow = 100, ncol = 2)
  
  result <- twinsGeneratorM(hap_ref, hap_seeds, count_obj, batch_size = 2, 
                            prior_window_size = 50, verbose = FALSE)
  expect_equal(nrow(result), 100)
})