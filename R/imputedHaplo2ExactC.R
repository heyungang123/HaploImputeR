#' @title Impute Haplotypes with Exact Allele Counts
#'
#' @description 
#' Generate haplotypes from a probability matrix while ensuring exact allele counts
#' match specified target values. Uses entropy-based ordering to determine the
#' imputation sequence for optimal matching.
#'
#' @details
#' The algorithm ensures exact allele count matching through:
#' \enumerate{
#'   \item Calculate allele probabilities from haplotype probability matrix
#'   \item Order chromosomes by entropy (most constrained first)
#'   \item Assign alleles sequentially with probability adjustment
#'   \item Track remaining counts and update probabilities
#' }
#'
#' This function is the core imputation engine used by \code{\link{haploSeeds}}
#' and the twinsGenerator family of functions.
#'
#' @param haplotype_Y_ref A matrix of reference haplotype types (sites × chromosomes).
#'   Used to identify unique haplotype patterns.
#' @param haplotype_Prob A matrix of probabilities (chromosomes × haplotype_types).
#'   Number of columns must match number of unique haplotypes in \code{haplotype_Y_ref}.
#' @param count_alleles A matrix of target allele counts (sites × 2).
#'   Column 1: count of allele 0, Column 2: count of allele 1.
#' @param seed Random seed for reproducibility. Default: NULL.
#' @param check_values Validate 0/1 values. Default: FALSE.
#'
#' @return A matrix of imputed haplotypes:
#'   \itemize{
#'     \item Rows: sites (same as \code{nrow(haplotype_Y_ref)})
#'     \item Columns: chromosomes (same as \code{nrow(haplotype_Prob)})
#'     \item Values: 0 or 1
#'   }
#'
#' @seealso
#' \code{\link{probComputationM}} for generating probability matrices\cr
#' \code{\link{haploSeeds}} which uses this function
#'
#' @importFrom stats runif
#' @export
#'
#' @examples
#' # Create example data
#' set.seed(42)
#' hap_Y <- matrix(sample(0:1, 40, replace = TRUE), nrow = 2, ncol = 20)
#'
#' # Determine number of unique haplotypes
#' unique_hap <- unique(as.data.frame(t(hap_Y)))
#' num_unique <- nrow(unique_hap)
#'
#' # Create probability matrix (10 chromosomes x unique haplotypes)
#' prob_mat <- matrix(runif(10 * num_unique), nrow = 10, ncol = num_unique)
#' prob_mat <- prob_mat / rowSums(prob_mat)
#'
#' # Define target counts
#' count <- matrix(c(5, 5, 5, 5), nrow = 2, ncol = 2, byrow = TRUE)
#'
#' # Impute haplotypes
#' imputed <- imputedHaplo2ExactC(hap_Y, prob_mat, count, seed = 42)
#'
#' # Verify dimensions
#' dim(imputed)  # 2 x 10
#'
#' # Verify allele counts
#' result_counts <- t(apply(imputed, 1, function(x) c(sum(x == 0), sum(x == 1))))
#' all(result_counts == count)  # TRUE
imputedHaplo2ExactC <- function(haplotype_Y_ref, haplotype_Prob, count_alleles,
                                seed = NULL, check_values = FALSE) {
  
  .setSeedIfProvided(seed)
  
  haplotype_Y_ref <- .validateInput(haplotype_Y_ref, "haplotype_Y_ref", check_values)
  haplotype_Prob <- .validateInput(haplotype_Prob, "haplotype_Prob", FALSE)
  count_alleles <- .validateInput(count_alleles, "count_alleles", FALSE)
  
  unique_haplotypes <- t(unique(t(haplotype_Y_ref)))
  unique_haplotypes <- .handleUniqueDim(unique_haplotypes)
  
  num_site <- nrow(unique_haplotypes)
  num_unique_haplotypes <- ncol(unique_haplotypes)
  
  unique_haplotypes <- .sortUniqueHaplotypes(unique_haplotypes)
  
  num_chr <- nrow(haplotype_Prob)
  hist_k <- matrix(0, nrow = num_chr, ncol = num_site)
  imputed_haplotype <- matrix(0, nrow = num_site, ncol = num_chr)
  
  haplotype_Prob <- as.matrix(haplotype_Prob)
  
  if (ncol(haplotype_Prob) != num_unique_haplotypes) {
    stop(paste0("haplotype_Prob columns (", ncol(haplotype_Prob), 
                ") must match number of unique haplotypes (", num_unique_haplotypes, ")"))
  }
  
  for (i in 1:num_site) {
    allele_Prob <- haplotype_Prob %*% as.numeric(unique_haplotypes[i, ])
    allele_Prob <- pmax(allele_Prob, .Machine$double.eps)
    allele_Prob <- pmin(allele_Prob, 1 - .Machine$double.eps)
    
    entropy_sites <- allele_Prob * .safeLog(allele_Prob) + 
                     (1 - allele_Prob) * .safeLog(1 - allele_Prob)
    temp_idx <- order(entropy_sites)
    
    for (j in 1:num_chr) {
      Qvalue <- sum(allele_Prob)
      Cvalue <- count_alleles[i, 2]
      
      if (Qvalue > Cvalue && Qvalue > 0) {
        k <- Cvalue / Qvalue
        hist_k[j, i] <- k
        prob_allele1 <- k * allele_Prob[temp_idx[j]]
      } else if (Qvalue < Cvalue) {
        denom <- num_chr + 1 - j - Qvalue
        if (denom > 0) {
          k <- (Cvalue - Qvalue) / denom
          hist_k[j, i] <- k
          prob_allele1 <- allele_Prob[temp_idx[j]] + k * (1 - allele_Prob[temp_idx[j]])
        } else {
          hist_k[j, i] <- 1
          prob_allele1 <- 1
        }
      } else {
        hist_k[j, i] <- 1
        prob_allele1 <- allele_Prob[temp_idx[j]]
      }
      
      prob_allele1 <- pmax(0, pmin(1, prob_allele1))
      
      temp_rand <- runif(1)
      temp_site_status <- as.integer(temp_rand <= prob_allele1)
      
      imputed_haplotype[i, temp_idx[j]] <- temp_site_status
      
      count_alleles[i, 1] <- count_alleles[i, 1] - (1 - temp_site_status)
      count_alleles[i, 2] <- count_alleles[i, 2] - temp_site_status
      
      tempPosi_unconsisitant_haplo <- which(unique_haplotypes[i, ] != temp_site_status)
      haplotype_Prob[temp_idx[j], tempPosi_unconsisitant_haplo] <- 0
      
      tempPosi_consistant_haplo <- which(unique_haplotypes[i, ] == temp_site_status)
      sum_prob <- sum(haplotype_Prob[temp_idx[j], tempPosi_consistant_haplo])
      if (sum_prob > 0) {
        haplotype_Prob[temp_idx[j], tempPosi_consistant_haplo] <- 
          haplotype_Prob[temp_idx[j], tempPosi_consistant_haplo] / sum_prob
      }
      
      allele_Prob[temp_idx[j]] <- 0
    }
  }
  
  return(imputed_haplotype)
}