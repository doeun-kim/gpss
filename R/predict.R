#' predict method for `gpss` objects
#'
#' @param object a model object for which prediction is desired.
#' @param newdata data frame on which to make predictions (test set)
#' @param type "response" or "scaled"
#' @param format "default" or "rvar"
#' @param interval "prediction" or "confidence"
#' @param level a numerical value between 0 and 1
#' @param prior_mean a numeric vector of prior mean values for Y at each test observation. Required when the model was trained with a \code{prior_mean}; see \code{\link{gp_predict}}. (default = NULL)
#' @param ... additional arguments (not used)
#' @inheritParams stats::predict
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
#' p <- predict(mod, newdata = dat_test)
#' p_confidence99 <- predict(mod, newdata = dat_test, interval = "confidence", level = 0.99)
#' @export
predict.gpss <- function(object, newdata = NULL, type = "response",
                         format = "default", interval = "confidence",
                         level = 0.95, prior_mean = NULL, ...) {

 # Input validation
  if (!isTRUE(format %in% c("default", "rvar"))) {
    stop('`format` must be "default" or "rvar".', call. = FALSE)
  }

  if (!isTRUE(interval %in% c("prediction", "confidence"))) {
    stop('`interval` must be "prediction" or "confidence".', call. = FALSE)
  }

  if (!isTRUE(level < 1 && level > 0)) {
    stop('`level` must be a numerical value between 0 and 1.', call. = FALSE)
  }

  type <- match.arg(type, c("response", "scaled"))

  if (!inherits(newdata, "data.frame")) {
    stop("`newdata` must be a data frame.", call. = FALSE)
  }

  if (!identical(type, "response") && identical(format, "rvar")) {
    stop('`format="rvar"` is only available with `type="response"`.', call. = FALSE)
  }

  # Get formula and create design matrix
  fo <- attr(object, "formula")
  if (!inherits(fo, "formula")) {
    stop("The `gpss` object must have been created using the formula interface of the `gpss()` function.",
         call. = FALSE)
  }

  # Drop the left-hand side
  fo <- as.formula(paste("~", deparse(formula(fo)[[3]])))
  X <- model.matrix(fo, newdata)

  # Check for missing columns
  miss <- setdiff(colnames(object$X), colnames(X))
  if (length(miss) > 0) {
    stop(sprintf("Missing columns in the design matrix: %s", paste(miss, collapse = ", ")),
         call. = FALSE)
  }
  X <- X[, colnames(object$X), drop = FALSE]

  # Get predictions
  out <- gp_predict(object, Xtest = X, prior_mean = prior_mean)

  # Select output type
  if (type == "response") {
    fit <- out$Ys_mean_orig
    gp_cov <- if (interval == "prediction") out$Ys_cov_orig else out$f_cov_orig
  } else {
    fit <- out$Ys_mean_scaled
    gp_cov <- if (interval == "prediction") out$Ys_cov_scaled else out$f_cov
  }

  # Format output
  if (format == "rvar") {
    if (!requireNamespace("posterior", quietly = TRUE)) {
      stop("Please install the `posterior` package.", call. = FALSE)
    }
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Please install the `MASS` package.", call. = FALSE)
    }

    rv <- MASS::mvrnorm(1e3, fit, gp_cov)
    out <- newdata
    out$rvar <- posterior::rvar(rv)
  } else {
    yvar <- diag(gp_cov)
    z <- qnorm(level + (1 - level) / 2)
    lwr <- fit - z * sqrt(yvar)
    upr <- fit + z * sqrt(yvar)
    out <- cbind(fit = fit, lwr = lwr, upr = upr)
  }

  return(out)
}
