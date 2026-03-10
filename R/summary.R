#' summary method for objects of class gpss
#'
#' @param object gpss object
#' @examples
#' library(gpss)
#' data(lalonde)
#' # categorical variables must be encoded as factors
#' dat <- transform(lalonde, race_ethnicity = factor(race_ethnicity))
#' # train and test sets
#' idx <- sample(seq_len(nrow(dat)), 500)
#' dat_train <- dat[idx, ]
#' dat_test <- dat[-idx, ]
#' # sample of data for speed
#' mod <- gpss(re78 ~ nsw + age + educ + race_ethnicity, data = dat_train)
#' summary(mod)
#' @keywords internal
#' @noRd
summary.gpss <- function(object, ...) {
  cat("GPSS Model Summary\n")
  cat("==================\n\n")

  # Basic info
  cat("Model Information\n")
  cat("-----------------\n")
  cat("Formula: "); print(attr(object, "formula"))
  cat("Observations:", nrow(object$Y), "\n")
  cat("Covariates:", ncol(object$X), "\n")
  cat("Scaled:", object$scale, "\n\n")

  # Kernel info
  cat("Kernel Settings\n")
  cat("---------------\n")
  cat("Type:", object$kernel_type, "\n")
  cat("Bandwidth (b):", object$b, "\n")
  cat("Noise variance (s2):", object$s2, "\n")

  if (!is.null(object$period_original)) {
    cat("Period (original):", object$period_original, "\n")
    cat("Period (scaled):", object$period_scaled, "\n")
    if (!is.null(object$time_col)) {
      cat("Time column:", object$time_col, "\n")
    } else {
      cat("Time column: first continuous column (default)\n")
    }
  }
  cat("\n")

  # Mixed data info
  cat("Data Composition\n")
  cat("----------------\n")
  cat("Mixed data:", object$mixed_data, "\n")
  if (isTRUE(object$mixed_data)) {
    cat("Categorical columns:", paste(object$cat_columns, collapse = ", "), "\n")
    cat("Processed categorical indices:", paste(object$cat_num_processed, collapse = ", "), "\n")
    cat("Processed continuous indices:", paste(object$cont_num_processed, collapse = ", "), "\n")
  }
  cat("\n")

  # Usage hints
  cat("Usage\n")
  cat("-----\n")
  cat("Fitted values: fit$post_mean_orig\n")
  cat("Standard errors: sqrt(diag(fit$post_cov_orig))\n")
  cat("Predictions: predict(fit, newdata)\n")

  invisible(object)
}
