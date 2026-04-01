# predict method for `gpss` objects

predict method for `gpss` objects

## Usage

``` r
# S3 method for class 'gpss'
predict(
  object,
  newdata = NULL,
  type = "response",
  format = "default",
  interval = "confidence",
  level = 0.95,
  prior_mean = NULL,
  ...
)
```

## Arguments

- object:

  a model object for which prediction is desired.

- newdata:

  data frame on which to make predictions (test set)

- type:

  "response" or "scaled"

- format:

  "default" or "rvar"

- interval:

  "prediction" or "confidence"

- level:

  a numerical value between 0 and 1

- prior_mean:

  a numeric vector of prior mean values for Y at each test observation.
  Required when the model was trained with a `prior_mean`; see
  [`gp_predict`](https://doeunkim.org/gpss/reference/gp_predict.md).
  (default = NULL)

- ...:

  additional arguments (not used)

## Value

If `format = "default"`, a numeric matrix with columns `fit`, `lwr`, and
`upr` giving the posterior mean and the lower and upper bounds of the
requested interval for each row of `newdata`. If `format = "rvar"`, a
copy of `newdata` with an additional column `rvar` containing a
posterior random variable (class `rvar` from the posterior package).

## Examples

``` r
# \donttest{
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
p <- predict(mod, newdata = dat_test)
p_confidence99 <- predict(mod, newdata = dat_test, interval = "confidence", level = 0.99)
# }
```
