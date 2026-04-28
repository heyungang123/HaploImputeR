#' @title Slice haplotype data with boundary protection
#' @description Internal function to safely slice haplotype matrices
#' @param haplotype Haplotype matrix
#' @param start_idx Start index
#' @param end_idx End index
#' @return Sliced matrix
#' @keywords internal
.sliceHaplotype <- function(haplotype, start_idx, end_idx) {
  start_idx <- max(1, start_idx)
  end_idx <- min(nrow(haplotype), end_idx)
  if (start_idx > end_idx) {
    return(matrix(nrow = 0, ncol = ncol(haplotype)))
  }
  haplotype[start_idx:end_idx, , drop = FALSE]
}

#' @title Validate input data
#' @description Check input matrix validity, enforce matrix type, and validate values
#' @param x Input object
#' @param name Object name for error message
#' @param check_values Whether to check for 0/1 values (default: FALSE)
#' @return Matrix representation of input
#' @keywords internal
.validateInput <- function(x, name = "input", check_values = FALSE) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop(paste0(name, " must be a matrix or data.frame"))
  }
  x <- as.matrix(x)
  if (nrow(x) == 0) {
    stop(paste0(name, " has zero rows"))
  }
  if (check_values) {
    unique_vals <- unique(as.numeric(x))
    unique_vals <- unique_vals[!is.na(unique_vals)]
    if (length(unique_vals) > 0 && !all(unique_vals %in% c(0, 1))) {
      stop(paste0(name, " must contain only 0 and 1 values (found: ", 
                  paste(unique_vals, collapse = ", "), ")"))
    }
  }
  x
}

#' @title Handle unique haplotype dimension issue
#' @description Convert unique result to matrix if dimension degenerates
#' @param unique_haplo Unique haplotype result
#' @keywords internal
.handleUniqueDim <- function(unique_haplo) {
  if (is.null(dim(unique_haplo))) {
    unique_haplo <- matrix(unique_haplo, nrow = 1, ncol = length(unique_haplo))
  }
  unique_haplo
}

#' @title Safe log computation
#' @description Compute log with protection against zero values
#' @param x Input values
#' @keywords internal
.safeLog <- function(x) {
  log(pmax(x, .Machine$double.eps))
}

#' @title Safe division
#' @description Division with protection against zero denominator
#' @param numerator Numerator
#' @param denominator Denominator
#' @keywords internal
.safeDivide <- function(numerator, denominator) {
  if (denominator == 0) return(0)
  numerator / denominator
}

#' @title Sort unique haplotypes
#' @description Sort unique haplotypes by site values
#' @param unique_haplo Unique haplotype matrix
#' @keywords internal
.sortUniqueHaplotypes <- function(unique_haplo) {
  num_site <- nrow(unique_haplo)
  if (num_site == 1) return(unique_haplo)
  for (i in 1:num_site) {
    unique_haplo <- unique_haplo[, order(unique_haplo[num_site + 1 - i, ])]
  }
  unique_haplo
}

#' @title Count unique haplotypes (vectorized)
#' @description Count occurrences of each unique haplotype pattern
#' @param haplotype_Y_train Training haplotype matrix
#' @param unique_haplotype Unique haplotype patterns
#' @return Vector of counts for each unique haplotype
#' @keywords internal
.countUniqueHaplotypes <- function(haplotype_Y_train, unique_haplotype) {
  num_site <- nrow(unique_haplotype)
  num_unique <- ncol(unique_haplotype)
  
  if (num_site == 1) {
    return(colSums(haplotype_Y_train == unique_haplotype))
  }
  
  matches <- sapply(1:num_unique, function(j) {
    colSums(haplotype_Y_train == unique_haplotype[, j]) == num_site
  })
  if (is.matrix(matches)) {
    return(colSums(matches))
  } else {
    return(as.numeric(matches) * ncol(haplotype_Y_train))
  }
}

#' @title Set random seed safely
#' @description Set random seed if provided, ensuring reproducibility
#' @param seed Random seed value (NULL means no seed set)
#' @keywords internal
.setSeedIfProvided <- function(seed) {
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1 || seed < 0) {
      stop("seed must be a single non-negative numeric value or NULL")
    }
    set.seed(seed)
  }
}

#' @title Check parallel package availability
#' @description Check if parallel processing is available and valid
#' @param parallel Logical, whether parallel processing is requested
#' @param n_workers Number of workers for parallel processing
#' @return Logical indicating if parallel processing is available
#' @keywords internal
.checkParallelAvailable <- function(parallel, n_workers) {
  if (!parallel) return(FALSE)
  if (!requireNamespace("parallel", quietly = TRUE)) {
    warning("Package 'parallel' not available, falling back to sequential processing")
    return(FALSE)
  }
  if (n_workers < 1) {
    warning("n_workers must be >= 1, using 1 worker")
    return(FALSE)
  }
  return(TRUE)
}