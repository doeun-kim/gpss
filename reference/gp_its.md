# Gaussian Process Interrupted Time Series

Gaussian Process Interrupted Time Series

## Usage

``` r
gp_its(
  y,
  dates,
  date_treat,
  covariates = NULL,
  kernel_type = "gaussian",
  b = NULL,
  s2 = 0.3,
  period = NULL,
  auto_period = FALSE,
  time_col = NULL,
  scale = TRUE,
  optimize = TRUE,
  interval_type = c("prediction", "confidence"),
  alpha = 0.05,
  mixed_data = FALSE,
  cat_columns = NULL,
  placebo_check = FALSE,
  placebo_periods = NULL,
  verbose = FALSE
)
```

## Arguments

- y:

  numeric vector or ts object of outcome values

- dates:

  Date vector of observation dates

- date_treat:

  Date of treatment/intervention

- covariates:

  optional matrix/data.frame of covariates

- kernel_type:

  kernel type for GP

- b:

  bandwidth (NULL for automatic)

- s2:

  noise variance (default 0.3)

- period:

  period for periodic kernels

- auto_period:

  logical; auto-detect period from pre-treatment data

- time_col:

  character or numeric; which column is the time variable

- scale:

  logical; scale covariates

- optimize:

  logical; optimize s2 via MLE

- interval_type:

  "prediction" or "confidence"

- alpha:

  significance level

- mixed_data:

  logical; whether covariates contain categorical variables

- cat_columns:

  character vector of categorical column names

- placebo_check:

  logical; run placebo checks

- placebo_periods:

  number of placebo periods (NULL for automatic)

- verbose:

  logical; print progress messages

## Value

An object of class `"gp_its"`, a list containing:

- y:

  the outcome vector

- dates:

  the date vector

- date_treat:

  the treatment date

- y0_hat:

  counterfactual predictions for post-treatment periods

- estimates:

  data frame with columns `tau_t`, `tau_cum`, `tau_avg`

- se:

  data frame of standard errors

- ci:

  data frame of confidence interval bounds

- gp_model:

  the fitted GP model from
  [`gp_train`](https://doeunkim.org/gpss/reference/gp_train.md)

- placebo_estimates:

  placebo check results (if requested)

## Examples

``` r
# \donttest{
# Simulated interrupted time series
set.seed(42)
n <- 50
dates <- seq(as.Date("2020-01-01"), by = "month", length.out = n)
y <- rnorm(n) + c(rep(0, 30), rep(2, 20))
res <- gp_its(y, dates, date_treat = as.Date("2022-07-01"))
print(res)
#> GP-ITS Model
#> ---------------------------------------- 
#> Observations: 50 (pre: 30, post: 20)
#> Treatment date: 2022-07-01
#> Kernel: gaussian
#> Interval type: prediction (95%)
#> 
#> Average treatment effect (final period):
#>   tau_avg = 1.918 (SE = 0.921)
#>   95% CI: [0.113, 3.723]
# }
```
