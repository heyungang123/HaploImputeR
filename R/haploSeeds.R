#' @title Generate Haplotype Seeds
#'
#' @description 
#' Generate initial haplotype seeds based on reference population frequency distribution.
#' Seeds are generated to match exact allele counts specified in the count matrix.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Identifies unique haplotype patterns in the training data
#'   \item Calculates frequency of each pattern
#'   \item Generates seeds that match target allele counts
#' }
#'
#' @param haplotype_Y_train A matrix of training haplotypes (sites × chromosomes).
#' @param count_alleles A matrix of target allele counts (sites × 2).
#' @param seed Random seed for reproducibility. Default: NULL.
#' @param check_values Validate 0/1 values. Default: FALSE.
#'
#' @return A matrix of generated haplotype seeds:
#'   \itemize{
#'     \item Rows: sites (same as \code{haplotype_Y_train})
#'     \item Columns: chromosomes (equal to \code{sum(count_alleles[1,])})
#'   }
#'
#' @seealso
#' \code{\link{imputedHaplo2ExactC}} for the underlying generation function\cr
#' \code{\link{twinsGenerator}} which uses this function internally
#'
#' @export
#'
#' @examples
#' # Create training data
#' set.seed(42)
#' hap_train <- matrix(sample(0:1, 200, replace = TRUE), nrow = 10, ncol = 20)
#'
#' # Define target counts (10 sites, 10 target chromosomes)
#' count <- matrix(c(5, 5), nrow = 10, ncol = 2)
#'
#' # Generate seeds
#' seeds <- haploSeeds(hap_train, count, seed = 42)
#'
#' # Verify dimensions
#' dim(seeds)  # 10 x 10
#'
#' # Verify allele counts
#' result_counts <- t(apply(seeds, 1, function(x) c(sum(x == 0), sum(x == 1))))
#' all(result_counts == count)  # TRUE
haploSeeds <- function(haplotype_Y_train, count_alleles, seed = NULL, 
                       check_values = FALSE) {
  
  .setSeedIfProvided(seed)
  
  haplotype_Y_train <- .validateInput(haplotype_Y_train, "haplotype_Y_train", check_values)
  count_alleles <- .validateInput(count_alleles, "count_alleles", FALSE)
  
  unique_haplotype <- t(unique(t(haplotype_Y_train)))
  unique_haplotype <- .handleUniqueDim(unique_haplotype)
  
  num_site <- nrow(unique_haplotype)
  num_unique_haplotypes <- ncol(unique_haplotype)
  
  unique_haplotype <- .sortUniqueHaplotypes(unique_haplotype)
  
  unique_haplotype_count <- .countUniqueHaplotypes(haplotype_Y_train, unique_haplotype)
  
  total_count <- sum(unique_haplotype_count)
  if (total_count == 0) {
    stop("No matching haplotypes found in training data")
  }
  
  unique_haplotype_freq <- unique_haplotype_count / total_count
  objectSize <- sum(count_alleles[1, ])
  
  transposed_matrix <- t(unique_haplotype_freq)
  temp_ProMatrix <- matrix(rep(transposed_matrix, objectSize),
                           nrow = nrow(transposed_matrix) * objectSize,
                           ncol = ncol(transposed_matrix))
  
  output_haplotype <- imputedHaplo2ExactC(haplotype_Y_train, temp_ProMatrix, 
                                          count_alleles, seed = seed,
                                          check_values = check_values)
  
  return(output_haplotype)
}