#' @title Calculate Linkage Disequilibrium (LD) R-squared
#'
#' @description 
#' Calculate pairwise LD r-squared value between two SNPs (single nucleotide polymorphisms).
#' Uses the standard formula:
#' \deqn{r^2 = \frac{(p_{AB} - p_A p_B)^2}{p_A(1-p_A)p_B(1-p_B)}}
#' where \eqn{p_A} and \eqn{p_B} are allele frequencies, and \eqn{p_{AB}} is the
#' frequency of the haplotype carrying both alleles.
#'
#' @details
#' The function computes the standard measure of linkage disequilibrium between
#' two biallelic loci. The r² value ranges from 0 (no LD) to 1 (complete LD).
#'
#' Edge cases handled:
#' \itemize{
#'   \item Returns NA when one or both loci are monomorphic (all same allele)
#'   \item Returns NA when denominator is zero
#' }
#'
#' @param LD_data A two-column matrix or data.frame where:
#'   \itemize{
#'     \item Each row represents a chromosome/haplotype
#'     \item Column 1: alleles at locus A (0 or 1)
#'     \item Column 2: alleles at locus B (0 or 1)
#'   }
#'
#' @return A numeric value:
#'   \itemize{
#'     \item R² value in range [0, 1] for valid calculations
#'     \item NA if calculation is undefined (monomorphic loci)
#'   }
#'
#' @seealso
#' \code{\link{winSizeR2}} which uses LD calculations\cr
#' \code{\link[genetics]{LD}} from the genetics package for more comprehensive LD measures
#'
#' @export
#'
#' @references
#' Lewontin, R. C. (1964). The interaction of selection and linkage. I. General 
#' considerations; heterotic models. Genetics, 49(1), 49-67.
#'
#' @examples
#' # Example 1: Perfect LD (r² = 1)
#' data1 <- matrix(c(0, 0, 0, 0, 1, 1, 1, 1), ncol = 2)
#' calculateLD(data1)  # 1
#'
#' # Example 2: No LD (r² = 0)
#' set.seed(42)
#' data2 <- matrix(sample(0:1, 200, replace = TRUE), ncol = 2)
#' calculateLD(data2)  # Close to 0
#'
#' # Example 3: Monomorphic locus (returns NA)
#' data3 <- matrix(c(0, 0, 0, 0), ncol = 2)
#' calculateLD(data3)  # NA
#'
#' # Example 4: Typical case
#' data4 <- matrix(c(0, 0, 0, 1, 1, 0, 1, 1), ncol = 2)
#' calculateLD(data4)  # 0.111...
calculateLD <- function(LD_data) {
  if (!is.matrix(LD_data) && !is.data.frame(LD_data)) {
    stop("LD_data must be a matrix or data.frame")
  }
  if (ncol(LD_data) != 2) {
    stop("LD_data must have exactly 2 columns")
  }
  if (nrow(LD_data) == 0) {
    stop("LD_data has zero rows")
  }
  
  LD_data <- as.matrix(LD_data)
  total_rows <- nrow(LD_data)
  
  AB_idx <- which(LD_data[, 1] == 0 & LD_data[, 2] == 0)
  pAB <- length(AB_idx) / total_rows
  
  A_idx <- which(LD_data[, 1] == 0)
  pA <- length(A_idx) / total_rows
  
  B_idx <- which(LD_data[, 2] == 0)
  pB <- length(B_idx) / total_rows
  
  denominator <- pA * (1 - pA) * pB * (1 - pB)
  
  if (denominator == 0 || is.na(denominator)) {
    return(NA)
  }
  
  rsquare <- (pAB - pA * pB)^2 / denominator
  
  return(rsquare)
}