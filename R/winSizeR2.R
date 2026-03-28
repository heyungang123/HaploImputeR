#' @title Compute Dynamic Window Size Based on R-squared
#'
#' @description 
#' Calculate adaptive window size for each site based on local LD (Linkage Disequilibrium)
#' structure. Sites with stronger correlations to upstream sites get larger window sizes,
#' while sites in regions of low LD get smaller windows.
#'
#' @details
#' The algorithm:
#' \enumerate{
#'   \item For each site, calculate R² with all upstream sites within \code{sizeMax}
#'   \item Count sites where R² exceeds \code{thresHold}
#'   \item Add \code{sizeBias} to the count
#'   \item Limit to maximum of site index (cannot exceed current position)
#' }
#'
#' This function supports parallel processing for large datasets through the
#' \code{parallel} and \code{n_workers} parameters.
#'
#' @param haplotype A matrix of haplotype data (sites × chromosomes).
#' @param thresHold R² threshold for counting correlated sites. Default: 0.05.
#'   Higher values result in smaller windows.
#' @param sizeMax Maximum number of upstream sites to consider. Default: 500.
#' @param sizeBias Base window size added to correlation count. Default: 200.
#' @param parallel Logical. Use parallel processing. Default: FALSE.
#' @param n_workers Number of parallel workers. Default: NULL (uses all cores).
#' @param seed Random seed for reproducibility. Default: NULL.
#' @param check_values Validate 0/1 values. Default: FALSE.
#'
#' @return A matrix of window sizes (sites × 1):
#'   \itemize{
#'     \item Row i contains the recommended window size for site i
#'     \item First row is always 0 (no upstream sites)
#'   }
#'
#' @seealso
#' \code{\link{twinsGeneratorRD}} which uses dynamic window sizes\cr
#' \code{\link{calculateLD}} for pairwise LD calculation
#'
#' @export
#'
#' @examples
#' # Create example haplotype data
#' set.seed(42)
#' hap <- matrix(sample(0:1, 1000, replace = TRUE), nrow = 100, ncol = 10)
#'
#' # Calculate dynamic window sizes
#' win_sizes <- winSizeR2(hap, thresHold = 0.05, sizeMax = 50, sizeBias = 20)
#'
#' # Check dimensions
#' dim(win_sizes)  # 100 x 1
#'
#' # Summary of window sizes
#' summary(win_sizes)
#'
#' \dontrun{
#' # Parallel processing for large datasets
#' large_hap <- matrix(sample(0:1, 100000, replace = TRUE), nrow = 1000, ncol = 100)
#' win_sizes_parallel <- winSizeR2(large_hap, parallel = TRUE, n_workers = 4)
#' }
winSizeR2 <- function(haplotype, thresHold = 0.05, sizeMax = 500, sizeBias = 200,
                      parallel = FALSE, n_workers = NULL, seed = NULL,
                      check_values = FALSE) {
  
  .setSeedIfProvided(seed)
  
  haplotype <- .validateInput(haplotype, "haplotype", check_values)
  
  if (!is.numeric(thresHold) || thresHold < 0 || thresHold > 1) {
    stop("thresHold must be a numeric value in [0, 1]")
  }
  if (!is.numeric(sizeMax) || sizeMax <= 0) {
    stop("sizeMax must be a positive numeric value")
  }
  if (!is.numeric(sizeBias) || sizeBias < 0) {
    stop("sizeBias must be a non-negative numeric value")
  }
  
  num_site <- nrow(haplotype)
  
  result <- matrix(0, nrow = num_site, ncol = 1)
  
  use_parallel <- .checkParallelAvailable(parallel, if(is.null(n_workers)) 1 else n_workers)
  
  if (use_parallel) {
    if (is.null(n_workers)) {
      n_workers <- parallel::detectCores()
    }
    n_workers <- min(n_workers, num_site - 1)
    
    site_indices <- 2:num_site
    
    cor_func <- function(i) {
      site_start <- max(1, i - sizeMax)
      site_end <- i - 1
      
      site_current <- haplotype[i, , drop = TRUE]
      site_before <- haplotype[site_start:site_end, , drop = FALSE]
      
      site_current <- as.numeric(site_current)
      
      if (nrow(site_before) == 0 || is.null(nrow(site_before))) {
        return(min(sizeBias, i))
      }
      
      cor_results <- apply(site_before, 1, function(x) {
        r <- cor(x, site_current)
        if (is.na(r)) return(0)
        return(r^2)
      })
      
      temp_count <- sum(cor_results >= thresHold, na.rm = TRUE)
      return(min(temp_count + sizeBias, i))
    }
    
    results_list <- parallel::mclapply(site_indices, cor_func, mc.cores = n_workers)
    
    for (k in seq_along(site_indices)) {
      result[site_indices[k], 1] <- results_list[[k]]
    }
  } else {
    for (i in 2:num_site) {
      site_start <- max(1, i - sizeMax)
      site_end <- i - 1
      
      site_current <- haplotype[i, , drop = TRUE]
      site_before <- haplotype[site_start:site_end, , drop = FALSE]
      
      site_current <- as.numeric(site_current)
      
      if (nrow(site_before) == 0 || is.null(nrow(site_before))) {
        result[i, 1] <- min(sizeBias, i)
        next
      }
      
      cor_results <- apply(site_before, 1, function(x) {
        r <- cor(x, site_current)
        if (is.na(r)) return(0)
        return(r^2)
      })
      
      temp_count <- sum(cor_results >= thresHold, na.rm = TRUE)
      result[i, 1] <- min(temp_count + sizeBias, i)
    }
  }
  
  return(result)
}