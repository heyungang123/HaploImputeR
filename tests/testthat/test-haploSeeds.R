test_that("haploSeeds validates input", {
  hap_train <- matrix(sample(0:1, 200, replace = TRUE), nrow = 10)
  count <- matrix(c(5, 5), nrow = 10, ncol = 2)
  
  expect_error(haploSeeds("not matrix", count), "matrix or data.frame")
  expect_error(haploSeeds(hap_train, "not matrix"), "matrix or data.frame")
})

test_that("haploSeeds validates 0/1 values when check_values=TRUE", {
  hap_invalid <- matrix(c(0, 1, 2, 0), nrow = 2)
  count <- matrix(c(1, 1), nrow = 2, ncol = 2)
  
  expect_error(haploSeeds(hap_invalid, count, check_values = TRUE), "only 0 and 1")
})

test_that("haploSeeds returns correct dimensions", {
  hap_train <- matrix(sample(0:1, 200, replace = TRUE), nrow = 10)
  count <- matrix(c(5, 5), nrow = 10, ncol = 2)
  
  result <- haploSeeds(hap_train, count)
  expect_equal(nrow(result), 10)
})

test_that("haploSeeds reproducibility with seed", {
  hap_train <- matrix(sample(0:1, 200, replace = TRUE), nrow = 10)
  count <- matrix(c(5, 5), nrow = 10, ncol = 2)
  
  result1 <- haploSeeds(hap_train, count, seed = 42)
  result2 <- haploSeeds(hap_train, count, seed = 42)
  expect_equal(result1, result2)
})