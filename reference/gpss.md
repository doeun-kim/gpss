# gpss: Gaussian Processes in Social Science

The Gaussian Process (GP) combines a highly flexible non-linear
regression approach with rigorous handling of uncertainty. A key feature
of this approach is that, while it can be used to produce a conditional
expectation function representing the mode of a posterior conditional
distribution, it does not “choose a particular model fit” and construct
uncertainty estimates conditional on putting full faith in that model.
This is valuable because the uncertainty estimates reflect the lesser
knowledge we have at locations near, at, or beyond the edge of the
observed data, where results from other approaches would become highly
model-dependent. We first offer an accessible explanation of GPs, and
provide an implementation more suitable to social science inference
problems, which reduces the number of user-chosen hyperparameters from
three to zero. We then illustrate the settings in which GPs can be most
valuable: those where conventional approaches have poor properties due
to model-dependency/extrapolation in data-sparse regions. Specifically,
we demonstrate the usefulness of GPs in contexts where (i) treated and
control models are needed by these groups have poor covariate overlap;
(ii) regression discontinuity, which depends on model estimates taken at
or just beyond the edge of their supporting data; and (iii) interrupted
time-series designs, where models are fitted prior to an event by
extrapolated after it.

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
  (dafault = TRUE)

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

- <http://doeunkim.org/gpss/>

- Report bugs at <https://github.com/doeun-kim/gpss/issues>

## Author

**Maintainer**: Chad Hazlett <chazlett@ucla.edu>

Authors:

- Soonhong Cho <tnsehdtm@gmail.com>

- Doeun Kim <doeun2@ucla.edu>
  ([ORCID](https://orcid.org/0000-0003-4789-6599))

## Examples

``` r
library(gpss)
data(lalonde)

# categorical variables must be encoded as factors
dat <- transform(lalonde, race_ethnicity = factor(race_ethnicity))

# train and test sets
idx <- sample(seq_len(nrow(dat)), 500)
dat_train <- dat[idx, ]
dat_test <- dat[-idx, ]

# Fit model
mod <- gpss(re78 ~ nsw + age + educ + race_ethnicity, data = dat_train)

# predictions in the test set
p <- predict(mod, dat_test)
length(p)
#> [1] 6525
head(p)
#>           fit       lwr      upr
#> [1,] 21633.55 19317.595 23949.51
#> [2,] 11967.40  6057.383 17877.43
#> [3,] 21674.79 19411.887 23937.69
#> [4,] 18540.16 15892.851 21187.47
#> [5,] 20811.21 16667.037 24955.39
#> [6,] 20521.15  5671.156 35371.15
```
