### Helper functions

#' one_hot
#'
#' to convert categorical variables into multiple binary variables with 1 and 0.
#'
#' @param data data frame containing only categorical variables
#
#' @return data frame containing expanded form of binary variables
#
#' @importFrom stats model.matrix contrasts
#' @keywords internal
#' @noRd
one_hot <- function(data) {
  onehot_data <- data.frame(lapply(data.frame(data),as.factor))
  onehot_data <- model.matrix(~ ., onehot_data,
                              contrasts.arg = lapply(onehot_data[,,drop = F],
                                                     contrasts,
                                                     contrasts=FALSE))
  onehot_data <- onehot_data[, -1]

  return(onehot_data)
}

#' mixed_data_processing
#'
#' to convert categorical variables in a data set into multiple binary variables with 1 and 0.
#'
#' @param X data set or matrix
#' @param cat_columns a numerical or character vector that indicates categorical variables
#' @param Xtest (optional) if there is testing data set separate from the training data, include the testing data
#
#' @return data frame containing expanded form of categorical variables
#' @keywords internal
#' @noRd
mixed_data_processing <- function(X, cat_columns = NULL, Xtest = NULL
) {

  # No categorical data case
  if (is.null(cat_columns)) { # If not, the data is all continuous. No processing is needed.
    if (!is.null(Xtest)) {
      X_test_processed <- as.matrix(Xtest)
    } else {
      X_test_processed <- NULL
    }
    # Return a list with the original data (as matrices) and column indices. All columns are marked as 'continuous'.
    return(list(X_train = as.matrix(X),
                X_test = X_test_processed,
                cat_num_processed = NULL,
                cont_num_processed = 1:ncol(X))) # All columns are continuous
  }

  # If cat_columns is provided, figure out which columns it refers to
  if (is.numeric(cat_columns)) {
    cat_num_idx <- cat_columns # if it's numeric, use it directly as indices
  } else if (is.character(cat_columns)) {
    cat_num_idx <- which(colnames(X) %in% cat_columns) # if it's character, find the matching column index
  }
  if (length(cat_num_idx) == 0) { stop("Specified cat_columns not found in X.") } # if no columns were matched, stop with an error

  n_train <- nrow(X)
  all_cols <- 1:ncol(X) # all column indices
  cont_num_idx <- setdiff(all_cols, cat_num_idx) # continuous column indices

  if (!is.null(Xtest)) { # the logic changes depending on whether we have a test set
    # CASE 1. We have a test set (Xtest is provided)
    # combine X and Xtest so the one-hot encoder sees all possible levels from both datasets at the same time
    Xtest_start <- n_train + 1
    Xtmp <- rbind(X, Xtest) # Combine rows
    Xtest_end <- nrow(Xtmp)

    # separate the combined data into categorical and continuous parts
    Xtmp_cat_df <- Xtmp[, cat_num_idx, drop = FALSE]
    allx_cat_matrix <- one_hot(Xtmp_cat_df) # one-hot coding: make dummy for each level
    Xtmp_cont_matrix <- Xtmp[, cont_num_idx, drop = FALSE]
    X_processed <- cbind(allx_cat_matrix, Xtmp_cont_matrix) # recombine new one-hot columns with original continuous columns
    X_processed <- apply(X_processed, 2, as.numeric) # ensure the final combined matrix is all numeric.

    # split the processed data back into train and test
    X_train_processed <- X_processed[1:n_train, , drop = FALSE]
    X_test_processed <- X_processed[Xtest_start:Xtest_end, , drop = FALSE]

    # new categorical columns are just 1 to N, where N is the number of new dummy variables created by one_hot()
    cat_num_processed <- 1:ncol(allx_cat_matrix)
    # new continuous columns come *after* the dummy variables
    cont_num_processed <- if(length(cont_num_idx) > 0) (ncol(allx_cat_matrix) + 1):ncol(X_processed) else NULL

  } else {
    # CASE 2 (simpler). We DO NOT have a test set (Xtest is NULL)
    # Separate X into categorical and continuous parts
    X_train_cat_df <- X[, cat_num_idx, drop = FALSE]
    X_train_cont_matrix <- X[, cont_num_idx, drop = FALSE]
    allx_cat_matrix <- one_hot(X_train_cat_df) # call one_hot() on the training categorical data only
    X_train_processed <- cbind(allx_cat_matrix, X_train_cont_matrix) # recombine one-hot columns and continuous columns
    X_train_processed <- apply(X_train_processed, 2, as.numeric) # make sure it's a numeric matrix.
    X_test_processed <- NULL # no test set was provided, so the processed version is NULL.

    # identify new column indices
    cat_num_processed <- 1:ncol(allx_cat_matrix)
    cont_num_processed <- if(length(cont_num_idx) > 0) (ncol(allx_cat_matrix) + 1):ncol(X_train_processed) else NULL
  }

  # this list contains everything the gp_train and gp_predict functions will need to handle the data correctly
  return(list(X_train = X_train_processed,     # the processed, one-hot encoded training data
              X_test = X_test_processed,      # the processed, one-hot encoded test data (or NULL)
              cat_num_processed = cat_num_processed, # indices of new dummy columns
              cont_num_processed = cont_num_processed # indices of new continuous columns
  ))
}

# getb_maxvar
#
#' @param X data set or matrix
#
#' @return cholesky-decomposed matrix
#
#' @keywords internal
#' @noRd
getb_maxvar <- function(x, kernel_type = "gaussian", period = NULL, maxsearch_b = 2000) {
  if (!is.numeric(maxsearch_b) || length(maxsearch_b) != 1) {
    stop("'maxsearch_b' must be a single numeric value.")
  }

  # Validate kernel_type matches available options
  kernel_type <- match.arg(kernel_type,
                           c("gaussian",
                             "gaussian_linear",
                             "gaussian_quadratic",
                             "gaussian_periodic_linear",
                             "gaussian_periodic_quadratic"))

  X <- as.matrix(x)
  n <- nrow(X)
  d <- ncol(X)

  # Validate period parameter for periodic kernels
  needs_period <- grepl("periodic", kernel_type)
  if (needs_period && (is.null(period) || !is.numeric(period) || period <= 0)) {
    stop(sprintf("Error: Period parameter required for kernel type '%s' and must be a positive numeric value.", kernel_type))
  }

  var_K <- function(b) {
    tryCatch({
      K <- kernel_symmetric(X, b, kernel_type, period)
      stats::var(K[lower.tri(K)])
    }, error = function(e) {
      return(0)  # Return 0 variance if computation fails
    })
  }
  res <- optimize(var_K, interval = c(0.01, maxsearch_b), maximum = TRUE)
  b_opt <- res$maximum

  return(b_opt)
}

#' kernel
#'
#' @param x1 first input matrix
#' @param x2 second input matrix
#' @param b bandwidth
#' @param kernel_type kernel type string
#' @param period period for periodic kernels
#' @return kernel matrix
#
#' @importFrom Rcpp sourceCpp
#' @keywords internal
#' @noRd
kernel <- function(x1, x2, b, kernel_type = "gaussian", period = NULL) {
  # Input validation
  if (missing(b) || is.null(b) || !is.numeric(b) || b <= 0) {
    stop("Error: 'b' must be a positive numeric value.")
  }

  kernel_type <- match.arg(kernel_type,
                           c("gaussian",
                             "gaussian_linear",
                             "gaussian_quadratic",
                             "gaussian_periodic_linear",
                             "gaussian_periodic_quadratic"))

  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)

  # Validate period parameter
  needs_period <- grepl("periodic", kernel_type)
  if (needs_period && (is.null(period) || !is.numeric(period) || period <= 0)) {
    stop(sprintf("Error: Period parameter required for kernel type '%s' and must be a positive numeric value.", kernel_type))
  }

  switch(kernel_type,
         gaussian = kernel_gaussian(x1, x2, b),
         gaussian_linear = kernel_gaussian_linear(x1, x2, b),
         gaussian_quadratic = kernel_gaussian_quadratic(x1, x2, b),
         gaussian_periodic_linear = kernel_gaussian_periodic_linear(x1, x2, b, period),
         gaussian_periodic_quadratic = kernel_gaussian_periodic_quadratic(x1, x2, b, period)
  )
}

#' kernel_symmetric
#'
#' @param x input matrix
#' @param b bandwidth
#' @param kernel_type kernel type string
#' @param period period for periodic kernels
#' @return symmetric kernel matrix
#
#' @importFrom Rcpp sourceCpp
#' @keywords internal
#' @noRd
kernel_symmetric <- function(x, b, kernel_type = "gaussian", period = NULL) {
  # Input validation
  if (missing(b) || is.null(b) || !is.numeric(b) || b <= 0) {
    stop("Error: 'b' must be a positive numeric value.")
  }

  kernel_type <- match.arg(kernel_type,
                           c("gaussian",
                             "gaussian_linear",
                             "gaussian_quadratic",
                             "gaussian_periodic_linear",
                             "gaussian_periodic_quadratic"))

  x <- as.matrix(x)

  # Validate period parameter
  needs_period <- grepl("periodic", kernel_type)
  if (needs_period && (is.null(period) || !is.numeric(period) || period <= 0)) {
    stop(sprintf("Error: Period parameter required for kernel type '%s' and must be a positive numeric value.", kernel_type))
  }

  switch(kernel_type,
         gaussian = kernel_symmetric_gaussian(x, b),
         gaussian_linear = kernel_symmetric_gaussian_linear(x, b),
         gaussian_quadratic = kernel_symmetric_gaussian_quadratic(x, b),
         gaussian_periodic_linear = kernel_symmetric_gaussian_periodic_linear(x, b, period),
         gaussian_periodic_quadratic = kernel_symmetric_gaussian_periodic_quadratic(x, b, period)
  )
}

#' estimate_period
#'
#' @param y numeric vector of observations
#' @param min_obs minimum number of observations required
#' @return estimated period (numeric)
#
#' @importFrom stats fft residuals lm var
#' @keywords internal
#' @noRd
estimate_period <- function(y, min_obs=7) {
  n <- length(y)
  if (n < min_obs) {
    warning(sprintf("Warning: Only %d observations (minimum %d recommended). Returning default period = 1.",
                    n, min_obs))
    return(1)
  }

  # Detrend the data
  tryCatch({
    time_index <- seq_along(y)
    trend_fit <- stats::lm(y ~ time_index)
    y_detrended <- stats::residuals(trend_fit)
  }, error = function(e) {
    warning("Warning: Could not detrend data. Using raw values.")
    y_detrended <- y - mean(y)
  })

  # Compute periodogram
  tryCatch({
    # Frequencies (excluding 0)
    n_freq <- floor(n / 2)
    frequencies <- seq(0, 0.5, length.out = n_freq + 1)[-1]
    fft_result <- fft(y_detrended) # Fast Fourier Transform
    periodogram <- Mod(fft_result[2:(n_freq + 1)])^2 / n # Power spectral density (periodogram); only use positive frequencies
    idx_max <- which.max(periodogram) # Find dominant frequency
    freq_dominant <- frequencies[idx_max]
    period_estimate <- 1 / freq_dominant # Convert frequency to period
    return(period_estimate)

  }, error = function(e) {
    warning(sprintf("Warning: Period estimation failed: %s. Returning default period = 1.", e$message))
    return(1)
  })
}
