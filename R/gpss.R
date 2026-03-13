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
#' \donttest{
#' # -- Basic fitting and prediction -----------------------------
#' library(gpss)
#' data(lalonde)
#' dat <- transform(lalonde, race_ethnicity = factor(race_ethnicity))
#'
#' idx <- sample(seq_len(nrow(dat)), 500)
#' mod <- gpss(re78 ~ nsw + age + educ + race_ethnicity, data = dat[idx, ])
#' p   <- predict(mod, dat[-idx, ])
#' head(p)
#'
#' # -- G-computation for the ATT (LaLonde data) ----------------
#' # The LaLonde dataset pairs 185 treated units from the NSW
#' # job-training program with 2490 PSID control units.  The
#' # experimental benchmark ATT is about $1794, but the naive
#' # difference in means is severely biased (~-$15k) due to
#' # poor covariate overlap.
#' #
#' # Strategy: fit a GP on the control group to learn E[Y(0)|X],
#' # predict counterfactual Y(0) for the treated, and compute
#' # ATT = mean(Y_obs - Y0_hat) with a Neyman-style SE.
#' # See Cho, Kim, and Hazlett (2026, Political Analysis).
#'
#' # Use gp_train / gp_predict for direct control over covariates
#' cat_vars <- c("black", "hisp", "married", "nodegr", "u74", "u75")
#' all_vars <- c("age", "educ", "re74", "re75",
#'               "black", "hisp", "married", "nodegr", "u74", "u75")
#'
#' X <- lalonde[, all_vars]
#' Y <- lalonde$re78
#' D <- lalonde$nsw
#'
#' # Fit GP on control group (optimize s2 via MLE)
#' gp_mod <- gp_train(X = X[D == 0, ], Y = Y[D == 0],
#'                    optimize = TRUE,
#'                    mixed_data = TRUE, cat_columns = cat_vars)
#'
#' # Predict E[Y(0) | X] for treated units
#' gp_pred <- gp_predict(gp_mod, Xtest = X[D == 1, ])
#'
#' Y1     <- Y[D == 1]                    # observed Y(1)
#' Y0_hat <- gp_pred$Ys_mean_orig         # predicted E[Y(0)|X]
#' tau_i  <- Y1 - Y0_hat                  # unit-level effects
#' att    <- mean(tau_i)                   # ATT point estimate
#'
#' # -- Neyman-style SE for the ATT --
#' # Two sources of uncertainty:
#' #
#' # 1. Sampling variance of unit-level effects: var(tau_i) / n1
#' #    This captures heterogeneity in treatment effects across
#' #    units, including variation in observed Y(1).
#' #
#' # 2. Posterior uncertainty in the counterfactuals: because Y(0)
#' #    predictions are correlated across units through the GP
#' #    posterior, we need the full covariance matrix V0 = f_cov_orig
#' #    (the posterior covariance of E[Y(0)|X] at treated X values).
#' #    Its contribution to Var(ATT) is (1/n1^2) * sum(V0),
#' #    i.e. 1'V0 1 / n1^2.
#'
#' n1     <- sum(D)
#' V0     <- gp_pred$f_cov_orig
#' se_att <- sqrt(var(tau_i) / n1 + sum(V0) / n1^2)
#'
#' cat(sprintf("GP ATT:  $%.0f (SE = $%.0f)\n", att, se_att))
#' cat(sprintf("95%% CI: [$%.0f, $%.0f]\n",
#'             att - 1.96 * se_att, att + 1.96 * se_att))
#' # GP estimate is close to the experimental benchmark (~$1794);
#' # the wide CI reflects genuine extrapolation uncertainty.
#' }
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
