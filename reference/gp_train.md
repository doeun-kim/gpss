# gp_train

to train GP model with training data set

## Usage

``` r
gp_train(
  X,
  Y,
  b = NULL,
  s2 = 0.3,
  optimize = FALSE,
  scale = TRUE,
  kernel_type = "gaussian",
  period = NULL,
  time_col = NULL,
  mixed_data = FALSE,
  cat_columns = NULL,
  Xtest = NULL,
  prior_mean = NULL
)
```

## Arguments

- X:

  a set of covariate data frame or matrix

- Y:

  Y vector (outcome variable)

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

- kernel_type:

  a character value indicating the kernel type (default = "gaussian")

- period:

  a numeric value for the period parameter, required for periodic
  kernels (default = NULL)

- time_col:

  a character or numeric value indicating which column is the time
  variable. If specified, this column will be moved to the first
  position for correct period scaling in periodic kernels. (default =
  NULL)

- mixed_data:

  a logical value to indicate whether the covariates contain a
  categorical/binary variable (default = FALSE)

- cat_columns:

  a character or a numerical vector indicating categorical variables.
  Must be character (not numeric) when time_col is specified. (default =
  NULL)

- Xtest:

  a data frame or a matrix of testing covariates. This is necessary when
  a non-overlapping categorical value exists between training and
  testing data sets. (default = NULL)

- prior_mean:

  a numeric vector of prior mean values for Y at each training
  observation. If provided, the GP is fitted to the residuals (Y -
  prior_mean), and the prior mean is added back to recover predictions
  on the original Y scale. Must be the same length as Y. (default =
  NULL)

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

## Examples

``` r
data(lalonde)
cat_vars <- c("race_ethnicity", "married")
all_vars <- c("age","educ","re74","re75","married", "race_ethnicity")

X <- lalonde[,all_vars]
Y <- lalonde$re78
D <- lalonde$nsw

X_train <- X[D==0,]
Y_train <- Y[D==0]
X_test <- X[D == 1,]
Y_test <- Y[D == 1]

gp_train.out <- gp_train(X = X_train, Y = Y_train,
optimize=TRUE, mixed_data = TRUE,
cat_columns = cat_vars)
gp_predict.out <- gp_predict(gp_train.out, X_test)
```
