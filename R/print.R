#' print method for objects of class gpss
#'
#' @export
print.gpss <- function(x) {
  cat("gpss object contains information, including...\n")
  cat("formula: "); print(attr(x, "formula"))
  cat("b (bandwidth):", x$b, "\n")
  cat("s2 (sigma squared):", x$s2, "\n")
  cat("mixed data (containing a categorical variable?):", x$mixed_data, "\n")
  if(x$mixed_data == TRUE){
    cat("categorical columns:", x$cat_columns,"\n")
    cat("categorical column numbers:", x$cat_num,"\n\n")
  }else{
    cat("\n")
  }
  cat("Other available information: posterior mean (scaled and original), posterior covariance (scaled and original), a kernel matrix of X, original values of Y and X, etc.", "\n\n")
  cat("For more detailed information, please use `summary(gpss object)`")
}


