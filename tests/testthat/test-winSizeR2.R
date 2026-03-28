test_that("winSizeR2 returns correct dimensions", {
  hap <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10)
  result <- winSizeR2(hap, 0.05, 5, 2)
  expect_equal(dim(result), c(10, 1))
})

test_that("winSizeR2 validates parameters", {
  hap <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10)
  expect_error(winSizeR2(hap, -1, 5, 2), "numeric value in \\[0, 1\\]")
  expect_error(winSizeR2(hap, 0.05, -1, 2), "positive numeric value")
  expect_error(winSizeR2(hap, 0.05, 5, -1), "non-negative numeric value")
})

test_that("winSizeR2 handles edge cases", {
  # Single row
  hap <- matrix(sample(0:1, 10, replace = TRUE), nrow = 1)
  result <- winSizeR2(hap, 0.05, 5, 2)
  expect_equal(dim(result), c(1, 1))
  expect_equal(result[1, 1], 0)
})

test_that("winSizeR2 respects sizeBias limit", {
  hap <- matrix(sample(0:1, 1000, replace = TRUE), nrow = 100)
  result <- winSizeR2(hap, 0.05, 500, 200)
  # Window size should not exceed site index
  for (i in 2:nrow(result)) {
    expect_true(result[i, 1] <= i)
  }
})

test_that("winSizeR2 reproducibility with seed", {
  hap <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10)
  result1 <- winSizeR2(hap, 0.05, 5, 2, seed = 42)
  result2 <- winSizeR2(hap, 0.05, 5, 2, seed = 42)
  expect_equal(result1, result2)
})

test_that("winSizeR2 parallel processing works", {
  skip_if_not(requireNamespace("parallel", quietly = TRUE))
  hap <- matrix(sample(0:1, 1000, replace = TRUE), nrow = 100)
  result_seq <- winSizeR2(hap, 0.05, 50, 20, parallel = FALSE)
  result_par <- winSizeR2(hap, 0.05, 50, 20, parallel = TRUE, n_workers = 2)
  expect_equal(dim(result_seq), dim(result_par))
})