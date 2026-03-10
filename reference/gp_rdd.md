# gp_rdd

This function implements a regression discontinuity (RD) design using
Gaussian Process (GP) regression on each side of the cutoff. Separate GP
models are estimated for observations below and above the threshold,
allowing for flexible and fully nonparametric functional forms. The
treatment effect is computed as the difference between the predicted
conditional means at the cutoff from the right and left limits. Standard
errors and confidence intervals are constructed using the posterior
predictive variance from the two independently fitted GPs.

## Usage

``` r
gp_rdd(
  X,
  Y,
  cut,
  alpha = 0.05,
  b = NULL,
  trim = FALSE,
  trim_k_value = 0.1,
  scale = TRUE
)
```

## Arguments

- X:

  forcing variable

- Y:

  Y vector (outcome variable)

- cut:

  cut point

- alpha:

  confidence level (default = 0.05)

- b:

  bandwidth (default = NULL)

- trim:

  a logical value indicating whether you want to do an automatic trim at
  a specific value of trim_k_value (default=FALSE)

- trim_k_value:

  a numerical value indicating the kernel value that you want to trim
  above (default = 0.1)

- scale:

  a logical value indicating whether you want to scale the covariates
  (default = TRUE)

## Value

- tau:

  an estimated treatment effect

- se:

  the standard error of tau

## Examples

``` r
n <- 100
tau <- 3
cut <- 0
x <- rnorm(n, 0, 1)
y <- rnorm(n, 0, 1) + tau*ifelse(x>cut, 1, 0)
gp_rdd.out <- gp_rdd(x, y, cut)
gp_rdd_plot(gp_rdd.out)
```
