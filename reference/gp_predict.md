# gp_predict

to predict outcome values at testing points by feeding the results
obtained from gp_train()

## Usage

``` r
gp_predict(gp, Xtest, prior_mean = NULL)
```

## Arguments

- gp:

  a list-form object obtained from gp_train()

- Xtest:

  a data frame or a matrix of testing data set

- prior_mean:

  a numeric vector of prior mean values for Y at each test point.
  Required when the model was trained with a `prior_mean`; added to
  `Ys_mean_orig` to recover predictions on the original Y scale. Must be
  the same length as `nrow(Xtest)`. (default = NULL)

## Value

- Xtest_scaled:

  testing data in a scaled form

- Xtest:

  the original testing data set

- Ys_mean_scaled:

  the predicted values of Y in a scaled form

- Ys_mean_orig:

  the predicted values of Y in the original scale

- Ys_cov_scaled:

  covariance of predicted Y in a scaled form

- Ys_cov_orig:

  covariance of predicted Y in the original scale

- f_cov_orig:

  covariance of target function in the original scale

- b:

  the bandwidth value obtained from gp_train()

- s2:

  the s2 value obtained from gp_train()
