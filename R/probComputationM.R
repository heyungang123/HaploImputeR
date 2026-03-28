#' @title Compute Probability Matrix via Lasso Regression
#'
#' @description 
#' Compute haplotype probability matrix using multinomial Lasso regression.
#' Uses reference haplotypes to train a model and predict probabilities for
#' target haplotypes based on prior site information.
#'
#' @details
#' The algorithm:
#' \enumerate{
#'   \item Extracts unique haplotype patterns from \code{haplotype_Y_ref}
#'   \item Trains a multinomial Lasso regression model using \code{haplotype_X_ref} as features
#'   \item Predicts haplotype probabilities for \code{haplotype_X_prior}
#' }
#'
#' This function is a core component of the haplotype generation pipeline,
#' providing probability estimates that guide the imputation process.
#'
#' @param haplotype_Y_ref A matrix of response haplotype types (sites × chromosomes).
#'   Each column represents a chromosome's alleles at the sites being predicted.
#' @param haplotype_X_ref A matrix of training features (prior_sites × chromosomes).
#'   Prior site information used to train the Lasso model.
#' @param haplotype_X_prior A matrix of prediction features (prior_sites × target_chromosomes).
#'   Prior information for which probabilities are predicted.
#' @param alpha Elastic net mixing parameter. Default: 1 (pure Lasso).
#'   Values in [0, 1]: 1 = Lasso, 0 = Ridge, intermediate = Elastic Net.
#' @param check_values Validate 0/1 values. Default: FALSE.
#'
#' @return A data frame of predicted haplotype probabilities:
#'   \itemize{
#'     \item Rows: target chromosomes (same as \code{ncol(haplotype_X_prior)})
#'     \item Columns: unique haplotype types
#'     \item Values: probabilities (row sums equal 1)
#'   }
#'
#' @import glmnet
#'
#' @seealso
#' \code{\link{imputedHaplo2ExactC}} which uses the probability matrix\cr
#' \code{\link[glmnet]{glmnet}} for the underlying Lasso implementation
#'
#' @export
#'
#' @examples
#' # Create example data
#' set.seed(42)
#' hap_Y <- matrix(sample(0:1, 200, replace = TRUE), nrow = 2, ncol = 100)
#' hap_X_ref <- matrix(sample(0:1, 500, replace = TRUE), nrow = 50, ncol = 100)
#' hap_X_prior <- matrix(sample(0:1, 100, replace = TRUE), nrow = 50, ncol = 20)
#'
#' # Compute probability matrix
#' prob_mat <- probComputationM(hap_Y, hap_X_ref, hap_X_prior)
#'
#' # Check dimensions
#' dim(prob_mat)  # 20 x (number of unique haplotypes)
#'
#' # Verify probabilities sum to 1
#' all(abs(rowSums(prob_mat) - 1) < 1e-10)  # TRUE
probComputationM <- function(haplotype_Y_ref, haplotype_X_ref, haplotype_X_prior,
                             alpha = 1, check_values = FALSE) {
  
  haplotype_Y_ref <- .validateInput(haplotype_Y_ref, "haplotype_Y_ref", check_values)
  haplotype_X_ref <- .validateInput(haplotype_X_ref, "haplotype_X_ref", check_values)
  haplotype_X_prior <- .validateInput(haplotype_X_prior, "haplotype_X_prior", check_values)
  
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required but not installed")
  }
  
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("alpha must be a numeric value in [0, 1]")
  }
  
  if (nrow(haplotype_X_ref) != nrow(haplotype_X_prior)) {
    stop("haplotype_X_ref and haplotype_X_prior must have the same number of rows")
  }
  
  unique_haplotypes <- unique(as.data.frame(t(haplotype_Y_ref)))
  unique_haplotypes <- t(as.matrix(unique_haplotypes))
  unique_haplotypes <- .handleUniqueDim(unique_haplotypes)
  
  num_site <- nrow(unique_haplotypes)
  num_unique_haplotypes <- ncol(unique_haplotypes)
  
  unique_haplotypes <- .sortUniqueHaplotypes(unique_haplotypes)
  
  size_haplotypes <- ncol(haplotype_Y_ref)
  Y_classes_ref <- matrix(0, nrow = num_unique_haplotypes, ncol = size_haplotypes)
  
  for (i in 1:size_haplotypes) {
    for (j in 1:num_unique_haplotypes) {
      if (all(haplotype_Y_ref[, i] == unique_haplotypes[, j])) {
        Y_classes_ref[j, i] <- 1
      }
    }
  }
  
  tryCatch({
    lasso_output <- glmnet::glmnet(t(haplotype_X_ref), t(Y_classes_ref), 
                                    family = "multinomial", alpha = alpha)
    min_lambda <- min(lasso_output$lambda)
    
    lasso_prediction <- predict(lasso_output, newx = t(haplotype_X_prior), 
                                type = "response", s = min_lambda)
    lasso_prediction <- as.data.frame(lasso_prediction)
    
    return(lasso_prediction)
  }, error = function(e) {
    stop(paste("glmnet fitting failed:", e$message))
  })
}