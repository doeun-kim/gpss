#' gp_rdd
#'
#' This function implements a regression discontinuity (RD) design
#' using Gaussian Process (GP) regression on each side of the cutoff.
#' Separate GP models are estimated for observations below and above the
#' threshold, allowing for flexible and fully nonparametric functional forms.
#' The treatment effect is computed as the difference between the predicted
#' conditional means at the cutoff from the right and left limits.
#' Standard errors and confidence intervals are constructed using the
#' posterior predictive variance from the two independently fitted GPs.
#'
#' @param X forcing variable
#' @param Y Y vector (outcome variable)
#' @param cut cut point
#' @param alpha confidence level (default = 0.05)
#' @param b bandwidth (default = NULL)
#' @param trim a logical value indicating whether you want to do an automatic trim at a specific value of trim_k_value (default=FALSE)
#' @param trim_k_value a numerical value indicating the kernel value that you want to trim above (default = 0.1)
#' @param scale a logical value indicating whether you want to scale the covariates (default = TRUE)
#' @examples
#' n <- 100
#' tau <- 3
#' cut <- 0
#' x <- rnorm(n, 0, 1)
#' y <- rnorm(n, 0, 1) + tau*ifelse(x>cut, 1, 0)
#' gp_rdd.out <- gp_rdd(x, y, cut)
#' gp_rdd_plot(gp_rdd.out)
#' @importFrom stats sd optimize complete.cases qnorm
#' @importFrom Rcpp sourceCpp
#' @return \item{tau}{an estimated treatment effect}
#' \item{se}{the standard error of tau}
#' @export
gp_rdd <- function(X, Y, cut, alpha=0.05, b=NULL,
                   trim=FALSE, trim_k_value=0.1, scale=TRUE){
  cutpoint <- c(cut, cut)

  b_left <- b
  b_right <- b

  na.ok <- complete.cases(X) & complete.cases(Y)
  X <- as.matrix(X[na.ok])
  Y <- as.numeric(Y[na.ok])

  X_left <- X[X<cut]
  X_right <- X[X>cut]
  Y_left <- Y[X<cut]
  Y_right <- Y[X>cut]

  #trimming: use only X with kernel value <0.1
  if(isTRUE(trim) & isTRUE(scale))
    stop("If `trim==TRUE`, scale must be FALSE")
  if(isTRUE(trim)){
    #1 set b automatically using maxvarK
    b_left <- getb_maxvar(X_left)
    b_right <- getb_maxvar(X_right)

    #2 trim the sample as to remove points farther than some X value
    #whose kernel value is `trim_k_value` (default: 0.1)
    trim_at_left <- cut - sqrt(-1*b_left*log(trim_k_value))
    trim_at_right <- cut + sqrt(-1*b_right*log(trim_k_value))

    Y_left <- Y_left[X_left>trim_at_left]
    Y_right <- Y_right[X_right<trim_at_right]
    X_left <- X_left[X_left>trim_at_left]
    X_right <- X_right[X_right<trim_at_right]
  }

  ## fit GP on left side of cutoff
  gp_train_l <- gp_train(X = X_left, Y = Y_left, b=b_left, scale = scale, optimize = TRUE)
  gp_pred_l <- gp_predict(gp_train_l, cutpoint)
  ## fit GP on right side of cutoff
  gp_train_r <- gp_train(X = X_right, Y = Y_right, b=b_right, scale = scale, optimize = TRUE)
  gp_pred_r <- gp_predict(gp_train_r, cutpoint)

  ## obtain estimate, se, and CI
  pred_l <- gp_pred_l$Ys_mean_orig[1]
  pred_r <- gp_pred_r$Ys_mean_orig[1]
  tau <- pred_r - pred_l
  se <- sqrt(diag(gp_pred_l$f_cov_orig)[1] + diag(gp_pred_r$f_cov_orig)[1])
  ci <- c(lower=tau - qnorm(1-alpha/2)*se, upper=tau + qnorm(1-alpha/2)*se)

  results <- list(tau=tau, se=se, ci=ci,
                  pred_l=pred_l, pred_r=pred_r,
                  b_left=gp_train_l$b, b_right=gp_train_r$b,
                  s2_left=gp_train_l$s2, s2_right=gp_train_r$s2,
                  n_left=length(X_left), n_right=length(X_right), trim=trim,
                  X = X,
                  Y = Y,
                  gp_train_l = gp_train_l,
                  gp_train_r = gp_train_r,
                  cut = cut,
                  alpha = alpha)

  if(trim==TRUE){
    results <- append(results,
                      c(trim_at_left=trim_at_left,
                        trim_at_right=trim_at_right))
  }
  return(results)
}

#' gp_rdd_plot
#'
#' to draw an RD plot using the results obtained from gp_rdd()
#'
#' @param gp_rdd_res a list-form results obtained from gp_rdd()
#' @param l_col a character value indicating the color of the left side of the cutoff point (default = "blue")
#' @param r_col a character value indicating the color of the right side of the cutoff point (default = "red")
#' @examples
#' library(ggplot2)
#' n <- 100
#' tau <- 3
#' cut <- 0
#' x <- rnorm(n, 0, 1)
#' y <- rnorm(n, 0, 1) + tau*ifelse(x>cut, 1, 0)
#' gp_rdd.out <- gp_rdd(x, y, cut)
#' gp_rdd_plot(gp_rdd.out) +
#'  geom_vline(xintercept = cut, lty="dashed")
#' @importFrom ggplot2 ggplot geom_point geom_line geom_ribbon theme_minimal aes
#' @importFrom rlang .data
#' @importFrom stats as.formula formula model.frame model.response
#' @return \item{gg}{an RD ggplot}
#' @export
gp_rdd_plot <- function(gp_rdd_res,
                        l_col = "blue",
                        r_col = "red"){

  ggdt <- data.frame(X = gp_rdd_res$X,
                     Y = gp_rdd_res$Y)

  Xtest_left <- seq(min(ggdt$X[ggdt$X<gp_rdd_res$cut]), gp_rdd_res$cut, length.out=100)
  gp_pred_l <- gp_predict(gp_rdd_res$gp_train_l, Xtest_left)
  l_se <- sqrt(diag(gp_pred_l$f_cov_orig))

  Xtest_right <- seq(gp_rdd_res$cut, max(ggdt$X[ggdt$X>gp_rdd_res$cut]), length.out=100)
  gp_pred_r <- gp_predict(gp_rdd_res$gp_train_r,
                          Xtest_right)
  r_se <- sqrt(diag(gp_pred_r$f_cov_orig))

  left_side <- data.frame(x = Xtest_left,
                          y = gp_pred_l$Ys_mean_orig,
                          l_lwr = gp_pred_l$Ys_mean_orig - qnorm(1-gp_rdd_res$alpha/2)*l_se,
                          l_upr = gp_pred_l$Ys_mean_orig + qnorm(1-gp_rdd_res$alpha/2)*l_se
  )

  right_side <- data.frame(x = Xtest_right,
                           y = gp_pred_r$Ys_mean_orig,
                           r_lwr = gp_pred_r$Ys_mean_orig - qnorm(1-gp_rdd_res$alpha/2)*r_se,
                           r_upr = gp_pred_r$Ys_mean_orig + qnorm(1-gp_rdd_res$alpha/2)*r_se
  )

  gg <- ggplot2::ggplot() +
    ggplot2::geom_point(data = ggdt,
                        ggplot2::aes(x = .data$X, y = .data$Y),
                        alpha = 0.5) +
    ggplot2::geom_line(data = left_side,
                       ggplot2::aes(x = .data$x, y = .data$y),
                       col = l_col) +
    ggplot2::geom_ribbon(data = left_side,
                         ggplot2::aes(x = .data$x, ymin = .data$l_lwr, ymax = .data$l_upr),
                         col = NA, alpha = 0.3, fill = l_col) +
    ggplot2::geom_line(data = right_side,
                       ggplot2::aes(x = .data$x, y = .data$y),
                       col = r_col) +
    ggplot2::geom_ribbon(data = right_side,
                         ggplot2::aes(x = .data$x, ymin = .data$r_lwr, ymax = .data$r_upr),
                         col = NA, alpha = 0.3, fill = r_col) +
    ggplot2::theme_minimal()

  return(gg)
}


#' Print method for gp_rdd objects
#' @param x a gp_rdd object
#' @param ... additional arguments (not used)
#' @keywords internal
#' @export
print.gp_rdd <- function(x, ...) {
  cat("GP-RDD Model\n")
  cat(strrep("-", 40), "\n")
  cat(sprintf("Cut point: %.2f\n", x$cut))
  cat(sprintf("Trimmed?: %s\n", x$trim))
  cat(sprintf("Tau: %.3f\n", x$tau))
  cat(sprintf("Confidence Interval: [%.3f, %.3f]\n", x$ci[1], x$ci[2]))
  cat(sprintf("No. of observations on the left: %.0f\n", x$n_left))
  cat(sprintf("No. of observations on the right: %.0f\n", x$n_right))

  invisible(x)
}

#' Summary method for gp_rdd objects
#' @param object a gp_rdd object
#' @param ... additional arguments (not used)
#' @keywords internal
#' @export
summary.gp_rdd <- function(object, ...) {
  cat("GP-RDD Model Summary\n")
  cat(strrep("=", 50), "\n\n")

  cat("DATA\n")
  cat(sprintf("  Total observations: %d\n", length(object$Y)))
  cat(sprintf("  Observations left of cutoff: %d\n", object$n_left))
  cat(sprintf("  Observations right of cutoff: %d\n", object$n_right))
  cat(sprintf("  Cut point: %.3f\n", object$cut))
  cat(sprintf("  Trimmed?: %s\n", object$trim))
  if (!is.null(object$trim_at_left)) {
    cat(sprintf("  Trim left at: %.3f\n", object$trim_at_left))
    cat(sprintf("  Trim right at: %.3f\n", object$trim_at_right))
  }

  cat("\nMODEL\n")
  cat(sprintf("  Bandwidth left (b_left): %.4f\n", object$b_left))
  cat(sprintf("  Bandwidth right (b_right): %.4f\n", object$b_right))
  cat(sprintf("  Noise variance left (s2_left): %.4f\n", object$s2_left))
  cat(sprintf("  Noise variance right (s2_right): %.4f\n", object$s2_right))

  cat("\nTREATMENT EFFECT\n")
  cat(sprintf("  Predicted left at cutoff: %.3f\n", object$pred_l))
  cat(sprintf("  Predicted right at cutoff: %.3f\n", object$pred_r))
  cat(sprintf("  Tau: %.3f\n", object$tau))
  cat(sprintf("  Standard error: %.3f\n", object$se))
  cat(sprintf("  Confidence Interval: [%.3f, %.3f]\n", object$ci[1], object$ci[2]))

  invisible(object)
}
