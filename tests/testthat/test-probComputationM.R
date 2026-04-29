test_that("probComputationM validates input", {
  hap_Y <- matrix(sample(0:1, 200, replace = TRUE), nrow = 2, ncol = 100)
  hap_X_ref <- matrix(sample(0:1, 5000, replace = TRUE), nrow = 50, ncol = 100)
  hap_X_prior <- matrix(sample(0:1, 1000, replace = TRUE), nrow = 50, ncol = 20)
  
  expect_error(probComputationM("not matrix", hap_X_ref, hap_X_prior), "matrix or data.frame")
})

test_that("probComputationM validates alpha parameter", {
  hap_Y <- matrix(sample(0:1, 200, replace = TRUE), nrow = 2, ncol = 100)
  hap_X_ref <- matrix(sample(0:1, 5000, replace = TRUE), nrow = 50, ncol = 100)
  hap_X_prior <- matrix(sample(0:1, 1000, replace = TRUE), nrow = 50, ncol = 20)
  
  expect_error(probComputationM(hap_Y, hap_X_ref, hap_X_prior, alpha = -1), "\\[0, 1\\]")
  expect_error(probComputationM(hap_Y, hap_X_ref, hap_X_prior, alpha = 2), "\\[0, 1\\]")
})

test_that("probComputationM returns data.frame", {
  hap_Y <- matrix(sample(0:1, 200, replace = TRUE), nrow = 2, ncol = 100)
  hap_X_ref <- matrix(sample(0:1, 5000, replace = TRUE), nrow = 50, ncol = 100)
  hap_X_prior <- matrix(sample(0:1, 1000, replace = TRUE), nrow = 50, ncol = 20)
  
  result <- probComputationM(hap_Y, hap_X_ref, hap_X_prior)
  expect_true(is.data.frame(result))
})

test_that("probComputationM validates dimension matching", {
  hap_Y <- matrix(sample(0:1, 200, replace = TRUE), nrow = 2, ncol = 100)
  hap_X_ref <- matrix(sample(0:1, 5000, replace = TRUE), nrow = 50, ncol = 100)
  hap_X_prior <- matrix(sample(0:1, 1500, replace = TRUE), nrow = 30, ncol = 50)
  
  expect_error(probComputationM(hap_Y, hap_X_ref, hap_X_prior), "same number of rows")
})