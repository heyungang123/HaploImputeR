#' @title Generate Twin Haplotypes from Reference Population
#'
#' @description 
#' Generate synthetic haplotypes based on reference population data while preserving 
#' allele frequency constraints. This function automatically generates initial haplotype 
#' seeds and extends them using Lasso-based probabilistic imputation.
#'
#' @details
#' The function works in three main phases:
#' \enumerate{
#'   \item Generate initial haplotype seeds using \code{\link{haploSeeds}}
#'   \item Iteratively extend haplotypes using a sliding window approach
#'   \item For each extension step, use Lasso regression to predict haplotype probabilities
#' }
#'
#' The algorithm ensures that the generated haplotypes match the exact allele counts
#' specified in \code{count_alleles_obj}.
#'
#' @param haplotype_ref A matrix of reference population haplotypes.
#'   \itemize{
#'     \item Rows represent sites (SNPs)
#'     \item Columns represent chromosomes
#'     \item Values must be 0 or 1
#'   }
#' @param count_alleles_obj A matrix specifying target population allele counts.
#'   \itemize{
#'     \item Rows represent sites (must match \code{haplotype_ref})
#'     \item Two columns: [count of allele 0, count of allele 1]
#'     \item Sum of each row determines number of target chromosomes
#'   }
#' @param seed_size Number of initial sites for automatic seed generation.
#'   Recommended value: 10-50 sites.
#' @param batch_size Number of sites to extend in each iteration.
#'   Recommended value: 2. Larger values may reduce accuracy.
#' @param prior_window_size Number of prior sites to consider for prediction.
#'   Larger values increase computational cost but may improve accuracy.
#' @param seed Random seed for reproducibility. Default is NULL (no seed set).
#'   Set a numeric value for reproducible results.
#' @param verbose Logical. Print progress messages during execution.
#'   Default is TRUE.
#' @param check_values Logical. Validate that input matrices contain only 0/1.
#'   Default is FALSE. Set TRUE for strict validation.
#'
#' @return A matrix of generated haplotypes with:
#'   \itemize{
#'     \item Rows: sites (same as \code{count_alleles_obj})
#'     \item Columns: chromosomes (equal to \code{sum(count_alleles_obj[1,])})
#'     \item Values: 0 or 1
#'   }
#'
#' @seealso
#' \code{\link{twinsGeneratorM}} for generation with external seeds\cr
#' \code{\link{twinsGeneratorRD}} for generation with dynamic window size\cr
#' \code{\link{haploSeeds}} for seed generation function
#'
#' @export
#'
#' @examples
#' # Create example reference population
#' set.seed(42)
#' hap_ref <- matrix(sample(0:1, 1000, replace = TRUE), nrow = 100, ncol = 10)
#'
#' # Define target population allele counts (100 sites, 10 target chromosomes)
#' count_obj <- matrix(c(5, 5), nrow = 100, ncol = 2)
#'
#' # Generate synthetic haplotypes
#' simu_hap <- twinsGenerator(
#'   haplotype_ref = hap_ref,
#'   count_alleles_obj = count_obj,
#'   seed_size = 10,
#'   batch_size = 2,
#'   prior_window_size = 50,
#'   seed = 42
#' )
#'
#' # Verify dimensions
#' dim(simu_hap)  # 100 x 10
#'
#' # Verify allele counts match
#' result_counts <- t(apply(simu_hap, 1, function(x) c(sum(x == 0), sum(x == 1))))
#' all(result_counts == count_obj)  # TRUE
twinsGenerator <- function(haplotype_ref, count_alleles_obj, seed_size, batch_size, 
                           prior_window_size, seed = NULL, verbose = TRUE, 
                           check_values = FALSE) {
  
  .setSeedIfProvided(seed)
  
  haplotype_ref <- .validateInput(haplotype_ref, "haplotype_ref", check_values)
  count_alleles_obj <- .validateInput(count_alleles_obj, "count_alleles_obj", FALSE)
  
  num_sites <- nrow(count_alleles_obj)
  num_chr_target <- sum(count_alleles_obj[1, ])
  
  if (seed_size > num_sites) {
    stop("seed_size cannot exceed number of sites in count_alleles_obj")
  }
  if (seed_size <= 0) {
    stop("seed_size must be positive")
  }
  if (batch_size <= 0 || batch_size >= num_sites) {
    stop("batch_size must be positive and less than num_sites")
  }
  if (prior_window_size <= 0) {
    stop("prior_window_size must be positive")
  }
  
  result <- matrix(0, nrow = num_sites, ncol = num_chr_target)
  
  seed_haplotypes <- haploSeeds(haplotype_ref[1:seed_size, ], 
                                count_alleles_obj[1:seed_size, ],
                                seed = seed, check_values = check_values)
  result[1:seed_size, ] <- seed_haplotypes
  currentSite_posi <- seed_size + 1
  
  if (verbose) message("Starting generation from site ", currentSite_posi)
  
  while (currentSite_posi <= (num_sites - batch_size)) {
    if (currentSite_posi <= prior_window_size) {
      haplo_prior_obj <- result[1:(currentSite_posi - 1), , drop = FALSE]
      haplo_prior_ref <- .sliceHaplotype(haplotype_ref, 1, currentSite_posi - 1)
      haplo_train_ref <- .sliceHaplotype(haplotype_ref, currentSite_posi, currentSite_posi + batch_size - 1)
    } else {
      start_idx <- currentSite_posi - prior_window_size
      haplo_prior_obj <- result[start_idx:(currentSite_posi - 1), , drop = FALSE]
      haplo_prior_ref <- .sliceHaplotype(haplotype_ref, currentSite_posi - prior_window_size, currentSite_posi - 1)
      haplo_train_ref <- .sliceHaplotype(haplotype_ref, currentSite_posi, currentSite_posi + batch_size - 1)
    }
    temp_Count <- count_alleles_obj[currentSite_posi:(currentSite_posi + batch_size - 1), , drop = FALSE]
    
    temp_ProbMatrix <- probComputationM(haplo_train_ref, haplo_prior_ref, haplo_prior_obj)
    temp_haplotypes <- imputedHaplo2ExactC(haplo_train_ref, temp_ProbMatrix, temp_Count,
                                           seed = seed, check_values = check_values)
    result[currentSite_posi:(currentSite_posi + batch_size - 1), ] <- temp_haplotypes
    currentSite_posi <- currentSite_posi + batch_size
    
    if (verbose) message("Extended to site ", currentSite_posi)
  }
  
  if (currentSite_posi <= num_sites) {
    start_idx <- max(1, currentSite_posi - prior_window_size)
    haplo_prior_obj <- result[start_idx:(currentSite_posi - 1), , drop = FALSE]
    haplo_prior_ref <- .sliceHaplotype(haplotype_ref, currentSite_posi - prior_window_size, currentSite_posi - 1)
    haplo_train_ref <- .sliceHaplotype(haplotype_ref, currentSite_posi, nrow(haplotype_ref))
    temp_Count <- count_alleles_obj[currentSite_posi:nrow(count_alleles_obj), , drop = FALSE]
    
    temp_ProbMatrix <- probComputationM(haplo_train_ref, haplo_prior_ref, haplo_prior_obj)
    temp_haplotypes <- imputedHaplo2ExactC(haplo_train_ref, temp_ProbMatrix, temp_Count,
                                           seed = seed, check_values = check_values)
    result[currentSite_posi:num_sites, ] <- temp_haplotypes
  }
  
  return(result)
}