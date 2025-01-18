#' predict method for `gpss` objects
#'
#' @export
predict.gpss <- function(object, newdata = NULL, type = "response") {
  type <- match.arg(type, c("response", "scaled"))
  if (!inherits(newdata, "data.frame")) {
    msg <- "`newdata` must be a data frame."
    stop(msg, call. = FALSE)
  }

  fo <- attr(object, "formula")
  if (!inherits(fo, "formula")) {
    msg <- "The `gpss` object must have been created using the formula interface of the `gpss()` function."
    stop(msg, call. = FALSE)
  }

  X <- model.matrix(fo, newdata)

  miss <- setdiff(colnames(object$X), colnames(X))
  if (length(miss) > 0) {
    msg <- sprintf("Missing columns in the design matrix: %s", paste(miss, collapse = ", "))
    stop(msg, call. = FALSE)
  }
  X <- X[, colnames(object$X), drop = FALSE]

  out <- gp_predict(object, Xtest = X)

  if (type == "response") {
    out <- out$Ys_mean_orig
  } else if (type == "scaled") {
    out <- out$Ys_mean_scaled
  }

  return(out)
}
