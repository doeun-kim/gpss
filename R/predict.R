#' predict method for `gpss` objects
#'
#' @param newdata data frame on which to make predictions (test set)
#' @param type "response" or "scaled"
#' @param format "default" or "rvar"
#' @inheritParams stats::predict
#'
#' @export
predict.gpss <- function(object, newdata = NULL, type = "response", format = "default") {
  if (!isTRUE(format %in% c("default", "rvar"))) {
    msg <- '`format` must be "default" or "rvar".'
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

  if (type == "response") {
    if (format == "rvar") {
      if (!requireNamespace("posterior")) {
        msg <- "Please install the `posterior` package."
        stop(msg, call. = FALSE)
      }
      # se <- sqrt(diag(gp_pred_l$f_cov_orig)[1] + diag(gp_pred_r$f_cov_orig)[1])
      mu <- out$Ys_mean_orig
      se <- sqrt(diag(out$f_cov_orig))
      rv <- lapply(seq_along(mu), function(i) rnorm(1e4, mu[i], se[i]))
      rv <- do.call(cbind, rv)
      rv <- posterior::rvar(rv)
      out <- newdata
      out$rvar <- posterior::rvar(rv)
    } else {
      out <- out$Ys_mean_orig
    }
  } else if (type == "scaled") {
    out <- out$Ys_mean_scaled
  }

  return(out)
}
