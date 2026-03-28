#' @title Generate Twin Haplotypes with Dynamic Window Size
#'
#' @description 
#' Generate synthetic haplotypes using a dynamic window size that adapts based on
#' local LD (Linkage Disequilibrium) structure. The window size varies per site
#' according to correlation patterns in the reference population.
#'
#' @details
#' This function uses site-specific window sizes, allowing the algorithm to:
#' \itemize{
#'   \item Use larger windows in regions of high LD
#'   \item Use smaller windows in regions of low LD
#'   \item Adapt to local haplotype structure
#' }
#'
#' Window sizes can be calculated using \code{\link{winSizeR2}}, which analyzes
#' R² correlations between sites.
#'
#' @param haplotype_ref A matrix of reference population haplotypes (sites × chromosomes).
#' @param haploSeeds_obj A matrix of initial haplotype seeds (seed_sites × target_chromosomes).
#' @param count_alleles_obj A matrix of target population allele counts (sites × 2).
#' @param batch_size Number of sites to extend in each iteration. Default: 2.
#' @param prior_window_size A numeric vector of window sizes with length equal to
#'   the number of sites. Use \code{\link{winSizeR2}} to generate this vector.
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
#' \code{\link{twinsGeneratorM}} for fixed window size\cr
#' \code{\link{winSizeR2}} for calculating dynamic window sizes
#'
#' @export
#'
#' @examples
#' # Create data
#' set.seed(42)
#' hap_ref <- matrix(sample(0:1, 5000, replace = TRUE), nrow = 50, ncol = 100)
#' seeds <- matrix(sample(0:1, 100, replace = TRUE), nrow = 5, ncol = 20)
#' counts <- matrix(c(10, 10), nrow = 50, ncol = 2)
#'
#' # Calculate dynamic window sizes
#' win_sizes <- winSizeR2(hap_ref, thresHold = 0.05, sizeMax = 30, sizeBias = 10)
#'
#' # Generate with dynamic windows
#' result <- twinsGeneratorRD(
#'   haplotype_ref = hap_ref,
#'   haploSeeds_obj = seeds,
#'   count_alleles_obj = counts,
#'   batch_size = 2,
#'   prior_window_size = win_sizes,
#'   seed = 42
#' )
#'
#' dim(result)  # 50 x 20
twinsGeneratorRD <- function(haplotype_ref, haploSeeds_obj, count_alleles_obj, 
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
  
  if (!is.numeric(prior_window_size) || length(prior_window_size) != num_sites) {
    stop("prior_window_size must be a numeric vector with length equal to num_sites")
  }
  if (any(prior_window_size < 0)) {
    stop("all elements in prior_window_size must be non-negative")
  }
  # Ensure minimum window size of 1 for valid processing
  prior_window_size <- pmax(prior_window_size, 1)
  if (seed_size > num_sites) {
    stop("haploSeeds_obj rows cannot exceed count_alleles_obj rows")
  }
  if (batch_size <= 0 || batch_size >= num_sites) {
    stop("batch_size must be positive and less than num_sites")
  }
  
  result <- matrix(0, nrow = num_sites, ncol = num_chr)
  result[1:seed_size, ] <- as.matrix(haploSeeds_obj)
  currentSite_posi <- seed_size + 1
  
  if (verbose) message("Starting generation from site ", currentSite_posi)
  
  while (currentSite_posi <= (num_sites - batch_size)) {
    current_window <- prior_window_size[currentSite_posi]
    
    if (currentSite_posi <= current_window) {
      haplo_prior_obj <- result[1:(currentSite_posi - 1), , drop = FALSE]
      haplo_prior_ref <- .sliceHaplotype(haplotype_ref, 1, currentSite_posi - 1)
      haplo_train_ref <- .sliceHaplotype(haplotype_ref, currentSite_posi, currentSite_posi + batch_size - 1)
    } else {
      start_idx <- currentSite_posi - current_window
      haplo_prior_obj <- result[start_idx:(currentSite_posi - 1), , drop = FALSE]
      haplo_prior_ref <- .sliceHaplotype(haplotype_ref, currentSite_posi - current_window, currentSite_posi - 1)
      haplo_train_ref <- .sliceHaplotype(haplotype_ref, currentSite_posi, currentSite_posi + batch_size - 1)
    }
    temp_Count <- count_alleles_obj[currentSite_posi:(currentSite_posi + batch_size - 1), , drop = FALSE]
    
    temp_ProbMatrix <- probComputationM(haplo_train_ref, haplo_prior_ref, haplo_prior_obj)
    temp_haplotypes <- imputedHaplo2ExactC(haplo_train_ref, temp_ProbMatrix, temp_Count,
                                           seed = seed, check_values = check_values)
    result[currentSite_posi:(currentSite_posi + batch_size - 1), ] <- temp_haplotypes
    currentSite_posi <- currentSite_posi + batch_size
    
    if (verbose) message("Extended to site ", currentSite_posi, " (window size: ", current_window, ")")
  }
  
  if (currentSite_posi <= num_sites) {
    current_window <- prior_window_size[currentSite_posi]
    start_idx <- max(1, currentSite_posi - current_window)
    haplo_prior_obj <- result[start_idx:(currentSite_posi - 1), , drop = FALSE]
    haplo_prior_ref <- .sliceHaplotype(haplotype_ref, currentSite_posi - current_window, currentSite_posi - 1)
    haplo_train_ref <- .sliceHaplotype(haplotype_ref, currentSite_posi, nrow(haplotype_ref))
    temp_Count <- count_alleles_obj[currentSite_posi:nrow(count_alleles_obj), , drop = FALSE]
    
    temp_ProbMatrix <- probComputationM(haplo_train_ref, haplo_prior_ref, haplo_prior_obj)
    temp_haplotypes <- imputedHaplo2ExactC(haplo_train_ref, temp_ProbMatrix, temp_Count,
                                           seed = seed, check_values = check_values)
    result[currentSite_posi:num_sites, ] <- temp_haplotypes
  }
  
  return(result)
}