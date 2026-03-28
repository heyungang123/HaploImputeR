test_that("calculateLD returns NA for degenerate cases", {
  expect_true(is.na(calculateLD(matrix(c(0,0,0,0), ncol=2))))
  expect_true(is.na(calculateLD(matrix(c(1,1,1,1), ncol=2))))
})

test_that("calculateLD computes correct r-squared", {
  data <- matrix(c(0,0,0,1,1,0,1,1), ncol=2)
  result <- calculateLD(data)
  expect_equal(result, 0.1111111, tolerance = 1e-6)
})

test_that("calculateLD handles mixed patterns", {
  data <- matrix(c(0,1,0,1), ncol=2)
  expect_equal(calculateLD(data), 1, tolerance = 1e-6)
})

test_that("calculateLD validates input", {
  expect_error(calculateLD(matrix(c(0,1), ncol=1)), "exactly 2 columns")
  expect_error(calculateLD(matrix(nrow=0, ncol=2)), "zero rows")
  expect_error(calculateLD("not a matrix"), "matrix or data.frame")
})