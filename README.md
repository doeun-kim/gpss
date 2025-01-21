
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gpss

<!-- badges: start -->

<!-- badges: end -->

The Gaussian Process (GP) combines a highly flexible non-linear
regression approach with rigorous handling of uncertainty. A key feature
of this approach is that, while it can be used to produce a conditional
expectation function representing the mode of a posterior conditional
distribution, it does not \`\`choose a particular model fit’’ and
construct uncertainty estimates conditional on putting full faith in
that model. This is valuable because the uncertainty estimates reflect
the lesser knowledge we have at locations near, at, or beyond the edge
of the observed data, where results from other approaches would become
highly model-dependent. We first offer an accessible explanation of GPs,
and provide an implementation more suitable to social science inference
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

## Installation

You can install the latest version by running:

``` r
devtools::install_github('doeun-kim/gpss')
```

### Troubleshooting installation

This version uses `Rcpp` extensively for speed reasons. These means you
need to have the right compilers on your machine.

#### Windows

If you are on Windows, you will need to install
[RTools](https://cran.r-project.org/bin/windows/Rtools/) if you haven’t
already. If you still are having difficulty with installing and it says
that the compilation failed, try installing it without support for
multiple architectures:

``` r
devtools::install_github('doeun-kim/gpss', args=c('--no-multiarch'))
```

#### Mac OSX

In order to compile the `C++` in this package, `RcppArmadillo` will
require you to have compilers installed on your machine. You may already
have these, but you can install them by running:

``` bash
xcode-select --install
```

You may need to install the `armadillo` software to compile the package. This can be done using [the Homebrew package manager](https://brew.sh/).

```bash
brew install armadillo
```

If you are having problems with this install on Mac OSX, specifically if
you are getting errors with either `lgfortran` or `lquadmath`, then try
open your Terminal and try the following:

``` bash
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

Also see section 2.16
[here](http://dirk.eddelbuettel.com/code/rcpp/Rcpp-FAQ.pdf)

## Examples

### Formula interface

``` r
library(gpss)
data(lalonde)
# categorical variables must be encoded as factors
dat <- transform(lalonde, race_ethnicity = factor(race_ethnicity))
# train and test sets
idx <- sample(seq_len(nrow(dat)), 500)
dat_train <- dat[idx, ]
dat_test <- dat[-idx, ]
# sample of data for speed
mod <- gpss(re78 ~ nsw + age + educ + race_ethnicity, data = dat_train)
# predictions in the test set
p <- predict(mod, dat_test)
length(p)
#> [1] 2175
head(p)
#> [1] 21656.54 20947.01 22602.46 15800.08 17971.55 10245.55
```

### Matrix interface

``` r
library(gpss)
data(lalonde)
cat_vars <- c("race_ethnicity", "married")
all_vars <-  c("age","educ","re74","re75","married", "race_ethnicity")

X <- lalonde[,all_vars]


Y <- lalonde$re78
D <- lalonde$nsw

X_train <- X[D==0,]
Y_train <- Y[D==0]

X_test <- X[D == 1,]
Y_test <- Y[D == 1]

gp_train.out <- gp_train(X = X_train, Y = Y_train, optimize=TRUE,
                         mixed_data = TRUE, 
                         cat_columns = cat_vars)

gp_predict.out <- gp_predict(gp_train.out, X_test)
```

## Optional: Accelerate R (faster matrix operations in R)

Mac users can see a significant speed up (5-10x) by using Apple’s native
Accelerate BLAS library (vecLib). Upgrade to the latest version of R and
RStudio, then follow the steps outlined
[here](https://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html#Which-BLAS-is-used-and-how-can-it-be-changed_003f):

``` bash
cd /Library/Frameworks/R.framework/Resources/lib

# for vecLib use
ln -sf libRblas.vecLib.dylib libRblas.dylib
```

For Window Users and to learn more details regarding installation and
reversion, please see [this blog
post](https://higgicd.github.io/posts/accelerating_r/).

Here is another Blog Post that explains about the benefits of using
Accelerate BLAS libraries:
[MacOS](https://mpopov.com/blog/2019/06/04/faster-matrix-math-in-r-on-macos/)
and [MacOS with
M1](https://mpopov.com/blog/2021/10/10/even-faster-matrix-math-in-r-on-macos-with-m1/).
