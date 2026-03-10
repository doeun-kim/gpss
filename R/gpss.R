formula_to_design <- function(formula, data) {
  X <- model.matrix(formula, data)

  # gp_train assumes there is no intercept/constant columns
  constants <- apply(X, 2, function(x) length(unique(x)) == 1)
  X <- X[, !constants, drop = FALSE]

  return(X)
}


#' Train a Gaussian Process Models
#'
#' @param formula an object of class "formula": a symbolic description of the model to be fitted
#' @param data a data frame containing the variables in the model.
#' @inheritParams gp_train
#' @inherit gp_train return
#'
#' @examples
#' library(gpss)
#' data(lalonde)
#'
#' # categorical variables must be encoded as factors
#' dat <- transform(lalonde, race_ethnicity = factor(race_ethnicity))
#'
#' # train and test sets
#' idx <- sample(seq_len(nrow(dat)), 500)
#' dat_train <- dat[idx, ]
#' dat_test <- dat[-idx, ]
#'
#' # Fit model
#' mod <- gpss(re78 ~ nsw + age + educ + race_ethnicity, data = dat_train)
#'
#' # predictions in the test set
#' p <- predict(mod, dat_test)
#' length(p)
#' head(p)
#' @param prior_mean a numeric vector of prior mean values for Y at each training observation. See \code{\link{gp_train}} for details. (default = NULL)
#' @export
gpss <- function(formula, data, b = NULL, s2 = 0.3, optimize = FALSE, scale = TRUE,
                 prior_mean = NULL) {
  # check user input and return an informative error if the user supplied incorrect objects
  sanity_formula_data(formula, data)

  if (!is.null(s2) && (!is.numeric(s2) || !isTRUE(length(s2) == 1))) {
    msg <- "`s2` must be `NULL` or a numeric of length 1."
    stop(msg, call. = FALSE)
  }

  if (!isTRUE(optimize) && !isFALSE(optimize)) {
    msg <- "`optimize` must be a logical flag."
  }

  if (!isTRUE(scale) && !isFALSE(scale)) {
    msg <- "`scale` must be a logical flag."
  }

  # fit model
  Y <- model.response(model.frame(formula, data))
  X <- formula_to_design(formula, data)
  out <- gp_train(X, Y, b = b, s2 = s2, optimize = optimize, scale = scale,
                  prior_mean = prior_mean)

  # output
  attr(out, "formula") <- formula
  class(out) <- "gpss"
  return(out)
}
