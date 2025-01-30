#' predict method for `gpss` objects
#'
#' @param newdata data frame on which to make predictions (test set)
#' @param type "response" or "scaled"
#' @param format "default" or "rvar"
#' @param interval "prediction" or "confidence"
#' @param level a numerical value between 0 and 1
#' @inheritParams stats::predict
#'
#' @export
predict.gpss <- function(object, newdata = NULL, type = "response", format = "default", interval = "prediction", level = 0.95) {
  if (!isTRUE(format %in% c("default", "rvar"))) {
    msg <- '`format` must be "default" or "rvar".'
    stop(msg, call. = FALSE)
  }

  if (!isTRUE(interval %in% c("prediction", "confidence"))) {
    msg <- '`interval` must be "prediction" or "confidence".'
    stop(msg, call. = FALSE)
  }

  if (!isTRUE(level < 1 & level>0)) {
    msg <- '`level` must be a numerical value between 0 and 1.'
    stop(msg, call. = FALSE)
  }

  type <- match.arg(type, c("response", "scaled"))
  if (!inherits(newdata, "data.frame")) {
    msg <- "`newdata` must be a data frame."
    stop(msg, call. = FALSE)
  }

  if (!identical(type, "response") && identical(format, "rvar")) {
    msg <- '`format="rvar"` is only available with `type="response"`.'
    stop(msg, call. = FALSE)
  }

  fo <- attr(object, "formula")
  if (!inherits(fo, "formula")) {
    msg <- "The `gpss` object must have been created using the formula interface of the `gpss()` function."
    stop(msg, call. = FALSE)
  }
  # Drop the left-hand side
  fo <- as.formula(paste("~", deparse(formula(fo)[[3]])))

  X <- model.matrix(fo, newdata)

  miss <- setdiff(colnames(object$X), colnames(X))
  if (length(miss) > 0) {
    msg <- sprintf("Missing columns in the design matrix: %s", paste(miss, collapse = ", "))
    stop(msg, call. = FALSE)
  }
  X <- X[, colnames(object$X), drop = FALSE]

  out <- gp_predict(object, Xtest = X)

  if (type == "response"){
    fit  <- out$Ys_mean_orig
    if (interval == "prediction"){
      gp_cov <- out$Ys_cov_orig
    } else if (interval == "confidence"){
      gp_cov <- out$f_cov_orig
    }
  } else if (type == "scaled"){
    fit <- out$Ys_mean_scaled
    if(interval == "prediction"){
      gp_cov <- out$Ys_cov_scaled
    } else if (interval == "confidence"){
      gp_cov <- out$f_cov
    }
  }

  if (format == "rvar") {
    if (!requireNamespace("posterior")) {
      msg <- "Please install the `posterior` package."
      stop(msg, call. = FALSE)
    }
    if (!requireNamespace("mvnfast")) {
      msg <- "Please install the `mvnfast` package."
      stop(msg, call. = FALSE)
    }

    # slow
    rv <- MASS::mvrnorm(1e3, fit, gp_cov)
    draws <- posterior::rvar(rv)

    out <- newdata
    out$rvar <- posterior::rvar(rv)
  } else if (format == "default"){
    yvar <- diag(gp_cov)
    level  <- level + (1-level)/2
    lwr <- fit - qnorm(level)*sqrt(yvar)
    upr <- fit + qnorm(level)*sqrt(yvar)
    out <- cbind(fit, lwr, upr)
  }

  return(out)
}
