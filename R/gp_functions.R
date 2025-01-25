### Main GP functions

#' gp_optimize
#'
#' to optimizes S2 using MLE. This function automatically runs when gp_train(..., optimize=TRUE).
#'
#' @param K kernel matrix (of covariates)
#' @param Y Y vector (outcome variable)
#'
#' @importFrom stats sd optimize
#' @importFrom Rcpp sourceCpp
#' @return \item{s2opt}{optimized s2 value}
#' \item{nlml}{the likelihood value at the optimal point}
#' @export
gp_optimize <- function(K, Y, optim.tol=0.1) {

  #define the objective function
  nlml <- function(K, Y, s2){ #b is fixed
    return(-1*log_marginal_likelihood_cpp(K, Y, s2))
  }

  #find the optimal s2 (by MLE)
  #since we always scale y, it is between 0 and 1; is 0.05 good for lower bound?
  opt <- optimize(nlml, interval=c(0.05, 1), maximum=FALSE, K=K, Y=Y, tol=optim.tol)
  results <- list(s2opt = opt$minimum, nlml = opt$objective)
  return(results)
}

#' gp_train
#'
#' to train GP model with training data set
#'
#' @param X a set of covariate data frame or matrix
#' @param Y Y vector (outcome variable)
#' @param b bandwidth (default = NULL)
#' @param s2 noise or a fraction of Y not explained by X (default = 0.3)
#' @param optimize a logical value to indicate whether an automatic optimized value of S2 should be used. If FALSE, users must define s2. (default = FALSE)
#' @param scale a logical value to indicate whether covariates should be scaled. (dafault = TRUE)
#' @param mixed_data a logical value to indicate whether the covariates contain a categorical/binary variable (default = FALSE)
#' @param cat_columns a character or a numerical vector indicating categorical variables (default = NULL)
#' @param Xtest a data frame or a matrix of testing covariates. This is necessary when a non-overlapping categorical value exists between training and testing data sets. (default = NULL)
#'
#' @return \item{post_mean_scaled}{posterior distribution of Y in a scaled form}
#' \item{post_mean_orig}{posterior distribution of Y in an original scale}
#' \item{post_cov_scaled}{posterior covariance matrix in a scaled form}
#' \item{post_cov_orig}{posterior covariance matrix in an original scale}
#' \item{K}{a kernel matrix of X}
#' \item{prior_mean_scaled}{prior distribution of mean in a scaled form}
#' \item{X.orig}{the original matrix or data set of X}
#' \item{X.init}{the original matrix or data set of X with categorical variables in an expanded form}
#' \item{X.init.mean}{the initial mean values of X}
#' \item{X.init.sd}{the initial standard deviation values of X}
#' \item{Y.init.mean}{the initial mean value of Y}
#' \item{Y.init.sd}{the initial standard deviation value of Y}
#' \item{K}{the kernel matrix of X}
#' \item{Y}{scaled Y}
#' \item{X}{scaled X}
#' \item{b}{bandwidth}
#' \item{s2}{sigma squared}
#' \item{alpha}{alpha value in Rasmussen and Williams (2006) p.19}
#' \item{L}{L value in Rasmussen and Williams (2006) p.19}
#' \item{mixed_data}{a logical value indicating whether X contains a categorical/binary variable}
#' \item{cat_columns}{a character or a numerical vector indicating the location of categorical/binary variables in X}
#' \item{cat_num}{a numerical vector indicating the location of categorical/binary variables in an expanded version of X}
#' \item{Xcolnames}{column names of X}
#' @examples
#' data(lalonde)
#' cat_vars <- c("race_ethnicity", "married")
#' all_vars <- c("age","educ","re74","re75","married", "race_ethnicity")
#'
#' X <- lalonde[,all_vars]
#' Y <- lalonde$re78
#' D <- lalonde$nsw
#'
#' X_train <- X[D==0,]
#' Y_train <- Y[D==0]
#' X_test <- X[D == 1,]
#' Y_test <- Y[D == 1]
#'
#' gp_train.out <- gp_train(X = X_train, Y = Y_train, optimize=TRUE, mixed_data = TRUE, cat_columns = cat_vars)
#' gp_predict.out <- gp_predict(gp_train.out, X_test)
#' @importFrom stats sd
#' @importFrom Rcpp sourceCpp
#'
#' @export
gp_train <- function(X, Y, b = NULL, s2 = 0.3, optimize = FALSE,
                     scale = TRUE, mixed_data = FALSE, cat_columns = NULL, Xtest = NULL){
  cat_num <- NULL

  ## pre-process outcome Y
  Y.init <- as.numeric(Y)
  Y.init.mean <- mean(Y.init)
  Y.init.sd <- sd(Y.init)
  Y <- scale(Y, center = Y.init.mean, scale = Y.init.sd)

  ## pre-process covariates X
  X <- as.matrix(X)
  X.orig  <- X
  n <- nrow(X)
  d <- ncol(X)
  if(mixed_data == TRUE){
    X_mix <- mixed_data_processing(X, cat_columns = cat_columns, Xtest=Xtest)
    X <- X_mix$X
    cat_num <- X_mix$cat_num
  }
  if(is.null(colnames(X))){
    colnames(X) <- paste("x", 1:d, sep = "")
  }
  X.init  <- X

  if (scale == TRUE & mixed_data == FALSE){
    X.init.mean <- apply(X.init, 2, mean)
    X.init.sd <- apply(X.init, 2, sd)
    if (sum(X.init.sd == 0) | sum(is.na(X.init.sd)) > 0) {
      stop("at least one column in X is a constant, please remove the constant(s)")
    }
    X <- scale(X, center = TRUE, scale = X.init.sd)
  } else if (scale == TRUE & mixed_data == TRUE){

    allx_cat <- sqrt(0.5)*X[, cat_num, drop= FALSE]
    X_cont <- X[, -cat_num, drop = FALSE]
    ### in this case, don't we have to keep both X.init.mean for categorical and continuous??
    X.init.mean = apply(X_cont, 2, mean)
    X.init.sd <- apply(X_cont, 2, sd)
    allx_cont <- scale(X_cont, center = TRUE, scale = X.init.sd)
    X <- as.matrix(cbind(allx_cat, allx_cont))
  } else {
    X.init.sd <- 1
    X.init.mean <- 0
  }

  ## b choice by maximizing variance of K when `b=null`
  if(is.null(b) & mixed_data == FALSE){
    #print("Choice of b left null; choosing b to maximize var(K)")
    b = getb_maxvar(X)
  } else if (is.null(b) & mixed_data == TRUE){
    # DK: check this part
    b = getb_maxvar(X)
    #b = 2*ncol(X)
  }

  ## Calculate GP
  K <- kernel_symmetric(X, b=b)

  # optimize s2
  if (isTRUE(optimize)) { #optimize s2, given b
    opt <- gp_optimize(K=K, Y=Y)
    s2 <- opt$s2opt
  } #otherwise, user-specified s2 is given (or default s2)

  ## new (directly use chol2inv to get K^-1): 1.33x faster
  L <- chol(K + diag(s2, nrow(K)))
  K_inv <- Matrix::chol2inv(L)
  m <- rep(0, nrow(X)) #zero mean prior
  a <- K_inv %*% (Y - m) #alpha with simplified computation using K^-1
  post_mean_scaled <- K %*% a
  post_cov_scaled <- K - K %*% K_inv %*% K
  prior_mean_scaled <- m


  #rescale Y to original
  post_mean_orig <- post_mean_scaled*Y.init.sd + Y.init.mean
  post_cov_orig <- Y.init.sd^2 * post_cov_scaled

  results <- list(post_mean_scaled = post_mean_scaled,
                  post_mean_orig = post_mean_orig,
                  post_cov_scaled = post_cov_scaled,
                  post_cov_orig = post_cov_orig,

                  prior_mean_scaled = prior_mean_scaled,
                  X.orig = X.orig,
                  X.init = X.init,
                  X.init.mean = X.init.mean,
                  X.init.sd = X.init.sd,
                  Y.init.mean = Y.init.mean,
                  Y.init.sd = Y.init.sd,
                  K = K,
                  Y = Y,
                  X = X,
                  b = b,
                  s2 = s2,
                  alpha=c(a),
                  L = L,
                  mixed_data = mixed_data,
                  cat_columns = cat_columns,
                  cat_num = cat_num,
                  Xcolnames = colnames(X))

  return(results)
}

#' gp_predict
#'
#' to predict outcome values at testing points by feeding the results obtained from gp_train()
#'
#' @param gp a list-form object obtained from gp_train()
#' @param Xtest a data frame or a matrix of testing data set
#' @examples
#' data(lalonde)
#' cat_vars <- c("race_ethnicity", "married")
#' all_vars <- c("age","educ","re74","re75","married", "race_ethnicity")
#'
#' X <- lalonde[,all_vars]
#' Y <- lalonde$re78
#' D <- lalonde$nsw
#'
#' X_train <- X[D==0,]
#' Y_train <- Y[D==0]
#' X_test <- X[D == 1,]
#' Y_test <- Y[D == 1]
#'
#' gp_train.out <- gp_train(X = X_train, Y = Y_train, optimize=TRUE, mixed_data = TRUE, cat_columns = cat_vars)
#' gp_predict.out <- gp_predict(gp_train.out, X_test)
#' @importFrom Rcpp sourceCpp
#' @return \item{Xtest_scaled}{testing data in a scaled form}
#' \item{Xtest}{the original testing data set}
#' \item{Ys_mean_scaled}{the predicted values of Y in a scaled form}
#' \item{Ys_mean_orig}{the predicted values of Y in the original scale}
#' \item{Ys_cov_scaled}{covariance of predicted Y in a scaled form}
#' \item{Ys_cov_orig}{covariance of predicted Y in the original scale}
#' \item{f_cov_orig}{covariance of target function in the original scale}
#' \item{b}{the bandwidth value obtained from gp_train()}
#' \item{s2}{the s2 value obtained from gp_train()}
#' @export
gp_predict <- function(gp, Xtest){
  mixed_data <- gp$mixed_data
  Xtest_init <- as.matrix(Xtest)

  if(!isTRUE(mixed_data)){
    Xtest <- sweep(Xtest_init, MARGIN=2, gp$X.init.mean, FUN = "-")
    Xtest <- sweep(Xtest, MARGIN=2, gp$X.init.sd, FUN = "/")
  } else {
    Xtest_mix <- mixed_data_processing(gp$X.orig,
                                       cat_columns = gp$cat_columns,
                                       Xtest=Xtest_init)
    Xtest <- Xtest_mix$X

    if(sum(gp$Xcolnames != colnames(Xtest))>0){
      print("Column names of training and testing data sets are not matching. It is probable that the categorical variables are causing issues.")
    }

    Xtest_cat <- sqrt(0.5)*Xtest[, gp$cat_num, drop= F]
    Xtest_cont <- Xtest[, -gp$cat_num, drop = F]
    Xtest_cont <- sweep(Xtest_cont, 2, gp$X.init.mean, FUN = "-")
    Xtest_cont <- sweep(Xtest_cont, 2, gp$X.init.sd, FUN = "/")
    Xtest <- as.matrix(cbind(Xtest_cat, Xtest_cont))
  }

  if (ncol(gp$X) != ncol(Xtest)) {
    print("ncol(Xtest) differs from ncol(Xtrain).")
    if(sum(!(colnames(Xtest) %in% gp$Xcolnames)) > 0){
      Xtest <- Xtest[,-c(!(colnames(Xtest) %in% gp$Xcolnames))]
    }else if(sum(!(gp$Xcolnames %in% colnames(Xtest)))){
      gp$X <- gp$X[,-c(!(gp$Xcolnames %in% colnames(Xtest)))]
    }
  }

  ## Calculate GP (Algorithm 2.1. in Rasmussen & Williams)
  Kss <- kernel(Xtest, Xtest, b=gp$b) #K_{**}
  Ks <- kernel(Xtest, gp$X, b=gp$b)

  Ys_mean_scaled <- c(Ks %*% gp$alpha)
  Ys_mean_orig = Ys_mean_scaled * gp$Y.init.sd + gp$Y.init.mean

  # we add s2*I to cov(f*) to compute cov(y*)
  f_cov <- solve(t(gp$L), t(Ks))
  f_cov <- Kss - crossprod(f_cov) #cov for target function
  Ys_cov_scaled <- f_cov + diag(gp$s2, nrow(f_cov)) #cov for new observation
  Ys_cov_orig <- gp$Y.init.sd^2 * Ys_cov_scaled #can be used for prediction interval
  f_cov_orig <- gp$Y.init.sd^2 * f_cov #can be used for confidence interval

  results <- list(Xtest_scaled = Xtest,
                  Xtest = Xtest_init,
                  #Ks = Ks, Kss = Kss,
                  Ys_mean_scaled = Ys_mean_scaled,
                  Ys_mean_orig = Ys_mean_orig,
                  Ys_cov_scaled = Ys_cov_scaled,
                  Ys_cov_orig = Ys_cov_orig,
                  f_cov_orig = f_cov_orig,
                  f_cov = f_cov,
                  b = gp$b,
                  s2 = gp$s2)
  return(results)
}

#' gp_rdd
#'
#' to perform RD using GP functions
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
                  cut = cut)

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
                          l_lwr = gp_pred_l$Ys_mean_orig - 1.96*l_se,
                          l_upr = gp_pred_l$Ys_mean_orig + 1.96*l_se
  )

  right_side <- data.frame(x = Xtest_right,
                           y = gp_pred_r$Ys_mean_orig,
                           r_lwr = gp_pred_r$Ys_mean_orig - 1.96*r_se,
                           r_upr = gp_pred_r$Ys_mean_orig + 1.96*r_se
  )

  gg <- ggplot2::ggplot() +
    ggplot2::geom_point(data = ggdt,
                        ggplot2::aes(X, Y), alpha = 0.5) +
    ggplot2::geom_line(data = left_side,
                       ggplot2::aes(x, y), col = l_col) +
    ggplot2::geom_ribbon(data = left_side,
                         ggplot2::aes(x, ymin = l_lwr, ymax = l_upr),
                         col = NA, alpha = 0.3, fill = l_col) +
    ggplot2::geom_line(data = right_side,
                       ggplot2::aes(x, y),
                       col = r_col) +
    ggplot2::geom_ribbon(data = right_side,
                         ggplot2::aes(x, ymin = r_lwr, ymax = r_upr),
                         col = NA, alpha = 0.3, fill = r_col) +
    ggplot2::theme_minimal()

  return(gg)
}
