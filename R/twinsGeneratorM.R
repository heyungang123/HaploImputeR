#' @title Generate Twin Haplotypes with External Seeds
#'
#' @description 
#' Generate synthetic haplotypes using externally provided haplotype seeds.
#' This function is useful when initial haplotype patterns need to be controlled
#' or when continuing from a previous generation run.
#'
#' @details
#' Unlike \code{\link{twinsGenerator}}, this function requires pre-generated seeds,
#' allowing more control over the initial haplotype structure. The seeds can be:
#' \itemize{
#'   \item Generated using \code{\link{haploSeeds}}
#'   \item Extracted from previous runs
#'   \item Custom-designed for specific scenarios
#' }
#'
#' The reference population (\code{haplotype_ref}) can have more chromosomes than
#' the target population. This design allows using a larger reference for model
#' training while generating a smaller target population.
#'
#' @param haplotype_ref A matrix of reference population haplotypes (sites × chromosomes).
#'   Can have more columns than \code{haploSeeds_obj}.
#' @param haploSeeds_obj A matrix of initial haplotype seeds (seed_sites × target_chromosomes).
#'   The number of columns determines the target population size.
#' @param count_alleles_obj A matrix of target population allele counts (sites × 2).
#'   Must have same number of rows as desired output.
#' @param batch_size Number of sites to extend in each iteration. Default: 2.
#' @param prior_window_size Number of prior sites for prediction window.
#' @param seed Random seed for reproducibility. Default: NULL.
#' @param verbose Print progress messages. Default: TRUE.
#' @param check_values Validate 0/1 values. Default: FALSE.
#'
#' @return A matrix of generated haplotypes:
#'   \itemize{
#'     \item Rows: sites (same as \code{count_alleles_obj})
#'     \item Columns: chromosomes (same as \code{haploSeeds_obj})
#'   }
#'
#' @seealso
#' \code{\link{twinsGenerator}} for automatic seed generation\cr
#' \code{\link{twinsGeneratorRD}} for dynamic window size
#'
#' @export
#'
#' @examples
#' # Create reference population (1000 chromosomes)
#' set.seed(42)
#' hap_ref <- matrix(sample(0:1, 5000, replace = TRUE), nrow = 50, ncol = 100)
#'
#' # Create seeds (20 chromosomes)
#' seeds <- matrix(sample(0:1, 100, replace = TRUE), nrow = 5, ncol = 20)
#'
#' # Define allele counts for 20 chromosomes
#' counts <- matrix(c(10, 10), nrow = 50, ncol = 2)
#'
#' # Generate with external seeds
#' result <- twinsGeneratorM(
#'   haplotype_ref = hap_ref,
#'   haploSeeds_obj = seeds,
#'   count_alleles_obj = counts,
#'   batch_size = 2,
#'   prior_window_size = 20,
#'   seed = 42
#' )
#'
#' dim(result)  # 50 x 20
twinsGeneratorM <- function(haplotype_ref, haploSeeds_obj, count_alleles_obj, 
                            batch_size, prior_window_size, seed = NULL,
                            verbose = TRUE, check_values = FALSE) {
  
  .setSeedIfProvided(seed)
  
  haplotype_ref <- .validateInput(haplotype_ref, "haplotype_ref", check_values)
  haploSeeds_obj <- .validateInput(haploSeeds_obj, "haploSeeds_obj", check_values)
  count_alleles_obj <- .validateInput(count_alleles_obj, "count_alleles_obj", FALSE)
  
  num_sites <- nrow(count_alleles_obj)
  num_chr <- ncol(haploSeeds_obj)
  seed_size <- nrow(haploSeeds_obj)
  
  # Note: haplotype_ref and haploSeeds_obj can have different number of columns
  # This design allows using a larger reference population for model training
  # while generating a smaller target population.
  # The result will have the same number of columns as haploSeeds_obj.
  
  if (seed_size > num_sites) {
    stop("haploSeeds_obj rows cannot exceed count_alleles_obj rows")
  }
  if (batch_size <= 0 || batch_size >= num_sites) {
    stop("batch_size must be positive and less than num_sites")
  }
  if (prior_window_size <= 0) {
    stop("prior_window_size must be positive")
  }
  
  result <- matrix(0, nrow = num_sites, ncol = num_chr)
  result[1:seed_size, ] <- as.matrix(haploSeeds_obj)
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