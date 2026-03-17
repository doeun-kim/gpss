# gpss: Gaussian Processes for Social Science

Provides Gaussian Process (GP) regression tools for social science
inference problems. GPs combine flexible nonparametric regression with
principled uncertainty quantification: rather than committing to a
single model fit, the posterior reflects lesser knowledge at the edge of
or beyond the observed data, where other approaches become highly
model-dependent. The package reduces user-chosen hyperparameters from
three to zero and supplies convenience functions for regression
discontinuity ('gp_rdd'), interrupted time-series ('gp_its'), and
general GP fitting ('gpss', 'gp_train', 'gp_predict'). Methods are
described in Cho, Kim, and Hazlett (2026)
[doi:10.1017/pan.2026.10032](https://doi.org/10.1017/pan.2026.10032) .

## Usage

``` r
gpss(
  formula,
  data,
  b = NULL,
  s2 = 0.3,
  optimize = FALSE,
  scale = TRUE,
  prior_mean = NULL
)
```

## Arguments

- formula:

  an object of class "formula": a symbolic description of the model to
  be fitted

- data:

  a data frame containing the variables in the model.

- b:

  bandwidth (default = NULL)

- s2:

  noise or a fraction of Y not explained by X (default = 0.3)

- optimize:

  a logical value to indicate whether an automatic optimized value of S2
  should be used. If FALSE, users must define s2. (default = FALSE)

- scale:

  a logical value to indicate whether covariates should be scaled.
  (default = TRUE)

- prior_mean:

  a numeric vector of prior mean values for Y at each training
  observation. See
  [`gp_train`](https://doeunkim.org/gpss/reference/gp_train.md) for
  details. (default = NULL)

## Value

- post_mean_scaled:

  posterior distribution of Y in a scaled form

- post_mean_orig:

  posterior distribution of Y in an original scale

- post_cov_scaled:

  posterior covariance matrix in a scaled form

- post_cov_orig:

  posterior covariance matrix in an original scale

- K:

  a kernel matrix of X

- prior_mean_scaled:

  prior distribution of mean in a scaled form

- X.orig:

  the original matrix or data set of X

- X.init:

  the original matrix or data set of X with categorical variables in an
  expanded form

- X.init.mean:

  the initial mean values of X

- X.init.sd:

  the initial standard deviation values of X

- Y.init.mean:

  the initial mean value of Y

- Y.init.sd:

  the initial standard deviation value of Y

- K:

  the kernel matrix of X

- Y:

  scaled Y

- X:

  scaled X

- b:

  bandwidth

- s2:

  sigma squared

- alpha:

  alpha value in Rasmussen and Williams (2006) p.19

- L:

  L value in Rasmussen and Williams (2006) p.19

- mixed_data:

  a logical value indicating whether X contains a categorical/binary
  variable

- cat_columns:

  a character or a numerical vector indicating the location of
  categorical/binary variables in X

- cat_num:

  a numerical vector indicating the location of categorical/binary
  variables in an expanded version of X

- time_col:

  the time column specification used (or NULL if not specified)

- Xcolnames:

  column names of X

- prior_mean:

  the prior mean vector supplied at training (or NULL)

## See also

Useful links:

- <https://doeunkim.org/gpss/>

- <https://github.com/doeun-kim/gpss>

- Report bugs at <https://github.com/doeun-kim/gpss/issues>

## Author

**Maintainer**: Chad Hazlett <chazlett@ucla.edu>

Authors:

- Soonhong Cho <tnsehdtm@gmail.com>

- Doeun Kim <doeun2@ucla.edu>
  ([ORCID](https://orcid.org/0000-0003-4789-6599))

## Examples

``` r
# \donttest{
# -- Basic fitting and prediction -----------------------------
library(gpss)
data(lalonde)
dat <- transform(lalonde, race_ethnicity = factor(race_ethnicity))

idx <- sample(seq_len(nrow(dat)), 500)
mod <- gpss(re78 ~ nsw + age + educ + race_ethnicity, data = dat[idx, ])
p   <- predict(mod, dat[-idx, ])
head(p)
#>             fit        lwr      upr
#> [1,] 19603.7203  17317.980 21889.46
#> [2,] 17517.0400  12644.866 22389.21
#> [3,] 21435.1890  19262.500 23607.88
#> [4,] 17832.7472  15209.686 20455.81
#> [5,] 13077.4249   8776.605 17378.24
#> [6,]   360.2888 -11517.823 12238.40

# -- G-computation for the ATT (LaLonde data) ----------------
# The LaLonde dataset pairs 185 treated units from the NSW
# job-training program with 2490 PSID control units.  The
# experimental benchmark ATT is about $1794, but the naive
# difference in means is severely biased (~-$15k) due to
# poor covariate overlap.
#
# Strategy: fit a GP on the control group to learn E[Y(0)|X],
# predict counterfactual Y(0) for the treated, and compute
# ATT = mean(Y_obs - Y0_hat) with a Neyman-style SE.
# See Cho, Kim, and Hazlett (2026, Political Analysis).

# Use gp_train / gp_predict for direct control over covariates
cat_vars <- c("black", "hisp", "married", "nodegr", "u74", "u75")
all_vars <- c("age", "educ", "re74", "re75",
              "black", "hisp", "married", "nodegr", "u74", "u75")

X <- lalonde[, all_vars]
Y <- lalonde$re78
D <- lalonde$nsw

# Fit GP on control group (optimize s2 via MLE)
gp_mod <- gp_train(X = X[D == 0, ], Y = Y[D == 0],
                   optimize = TRUE,
                   mixed_data = TRUE, cat_columns = cat_vars)

# Predict E[Y(0) | X] for treated units
gp_pred <- gp_predict(gp_mod, Xtest = X[D == 1, ])

Y1     <- Y[D == 1]                    # observed Y(1)
Y0_hat <- gp_pred$Ys_mean_orig         # predicted E[Y(0)|X]
tau_i  <- Y1 - Y0_hat                  # unit-level effects
att    <- mean(tau_i)                   # ATT point estimate

# -- Neyman-style SE for the ATT --
# Two sources of uncertainty:
#
# 1. Sampling variance of unit-level effects: var(tau_i) / n1
#    This captures heterogeneity in treatment effects across
#    units, including variation in observed Y(1).
#
# 2. Posterior uncertainty in the counterfactuals: because Y(0)
#    predictions are correlated across units through the GP
#    posterior, we need the full covariance matrix V0 = f_cov_orig
#    (the posterior covariance of E[Y(0)|X] at treated X values).
#    Its contribution to Var(ATT) is (1/n1^2) * sum(V0),
#    i.e. 1'V0 1 / n1^2.

n1     <- sum(D)
V0     <- gp_pred$f_cov_orig
se_att <- sqrt(var(tau_i) / n1 + sum(V0) / n1^2)

cat(sprintf("GP ATT:  $%.0f (SE = $%.0f)\n", att, se_att))
#> GP ATT:  $1907 (SE = $1760)
cat(sprintf("95%% CI: [$%.0f, $%.0f]\n",
            att - 1.96 * se_att, att + 1.96 * se_att))
#> 95% CI: [$-1542, $5357]
# GP estimate is close to the experimental benchmark (~$1794);
# the wide CI reflects genuine extrapolation uncertainty.
# }
```
