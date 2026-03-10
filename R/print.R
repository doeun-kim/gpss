#' print method for objects of class gpss
#'
#' @param x gpss object
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
#' print(mod)
#' @keywords internal
#' @noRd
print.gpss <- function(x, ...) {
  cat("gpss object\n")
  cat("===========\n")
  cat("Formula: "); print(attr(x, "formula"))
  cat("Kernel type:", x$kernel_type, "\n")
  cat("Bandwidth (b):", x$b, "\n")
  cat("Noise variance (s2):", x$s2, "\n")

  # Periodic kernel info
  if (!is.null(x$period_original)) {
    cat("Period (original):", x$period_original, "\n")
    cat("Period (scaled):", x$period_scaled, "\n")
    if (!is.null(x$time_col)) {
      cat("Time column:", x$time_col, "\n")
    }
  }

  # Mixed data info
  cat("Mixed data:", x$mixed_data, "\n")
  if (isTRUE(x$mixed_data)) {
    cat("Categorical columns:", x$cat_columns, "\n")
  }

  cat("\nUse `summary()` for more details.\n")
  invisible(x)
}
