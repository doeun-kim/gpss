#' summary method for objects of class gpss
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
#' summary(mod)
#' @export
summary.gpss <- function(x) {
  cat("Basic Model Information\n")
  cat("formula: "); print(attr(x, "formula"))
  cat("number of observations: ", dim(x$Y)[1], "\n")
  cat("number of covariates: ", dim(x$X)[2], "\n")
  cat("mixed data (containing a categorical variable?):", x$mixed_data, "\n")
  if(x$mixed_data == TRUE){
    cat("categorical columns:", x$cat_columns,"\n")
  }else{
    cat("\n")
  }

  cat("Hyperparrameters\n")
  cat("b (bandwidth):", x$b, "\n")
  cat("s2 (noise variance):", x$s2, "\n\n")

  cat("Scaling information\n")
  cat("scaled:", x$mixed_data, "\n\n")

  cat("Usage Example\n")
  cat("e.g. fit <- gpss(Y~X)")
  cat("to extract SEs of fitted values: sqrt(diag(fit$post_cov_orig))")
}


