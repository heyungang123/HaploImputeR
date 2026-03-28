test_that("imputedHaplo2ExactC validates input", {
  set.seed(123)
  hap_Y <- matrix(sample(0:1, 40, replace = TRUE), nrow = 4)
  unique_hap <- unique(as.data.frame(t(hap_Y)))
  num_unique <- nrow(unique_hap)
  
  prob_mat <- matrix(runif(10 * num_unique), nrow = 10, ncol = num_unique)
  prob_mat <- prob_mat / rowSums(prob_mat)
  count <- matrix(c(5, 5), nrow = 4, ncol = 2)
  
  expect_error(imputedHaplo2ExactC("not matrix", prob_mat, count), "matrix or data.frame")
})

test_that("imputedHaplo2ExactC returns correct dimensions", {
  set.seed(123)
  hap_Y <- matrix(sample(0:1, 40, replace = TRUE), nrow = 4)
  unique_hap <- unique(as.data.frame(t(hap_Y)))
  num_unique <- nrow(unique_hap)
  
  prob_mat <- matrix(runif(10 * num_unique), nrow = 10, ncol = num_unique)
  prob_mat <- prob_mat / rowSums(prob_mat)
  count <- matrix(c(5, 5), nrow = 4, ncol = 2)
  
  result <- imputedHaplo2ExactC(hap_Y, prob_mat, count)
  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 10)
})

test_that("imputedHaplo2ExactC reproducibility with seed", {
  set.seed(123)
  hap_Y <- matrix(sample(0:1, 40, replace = TRUE), nrow = 4)
  unique_hap <- unique(as.data.frame(t(hap_Y)))
  num_unique <- nrow(unique_hap)
  
  prob_mat <- matrix(runif(10 * num_unique), nrow = 10, ncol = num_unique)
  prob_mat <- prob_mat / rowSums(prob_mat)
  count <- matrix(c(5, 5), nrow = 4, ncol = 2)
  
  result1 <- imputedHaplo2ExactC(hap_Y, prob_mat, count, seed = 42)
  result2 <- imputedHaplo2ExactC(hap_Y, prob_mat, count, seed = 42)
  expect_equal(result1, result2)
})

test_that("imputedHaplo2ExactC contains only 0/1 values", {
  set.seed(123)
  hap_Y <- matrix(sample(0:1, 40, replace = TRUE), nrow = 4)
  unique_hap <- unique(as.data.frame(t(hap_Y)))
  num_unique <- nrow(unique_hap)
  
  prob_mat <- matrix(runif(10 * num_unique), nrow = 10, ncol = num_unique)
  prob_mat <- prob_mat / rowSums(prob_mat)
  count <- matrix(c(5, 5), nrow = 4, ncol = 2)
  
  result <- imputedHaplo2ExactC(hap_Y, prob_mat, count)
  unique_vals <- unique(as.numeric(result))
  expect_true(all(unique_vals %in% c(0, 1)))
})

test_that("imputedHaplo2ExactC validates dimension matching", {
  set.seed(123)
  hap_Y <- matrix(sample(0:1, 40, replace = TRUE), nrow = 4)
  
  # Wrong number of columns in prob_mat
  prob_mat <- matrix(runif(30), nrow = 10, ncol = 3)
  prob_mat <- prob_mat / rowSums(prob_mat)
  count <- matrix(c(5, 5), nrow = 4, ncol = 2)
  
  expect_error(imputedHaplo2ExactC(hap_Y, prob_mat, count), "must match number of unique haplotypes")
})