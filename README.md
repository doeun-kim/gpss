
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gpss

<!-- badges: start -->
<!-- badges: end -->

The goal of gpss is to …

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

If you are having problems with this install on Mac OSX, specifically if
you are getting errors with either `lgfortran` or `lquadmath`, then try
open your Terminal and try the following:

``` bash
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

Also see section 2.16
[here](http://dirk.eddelbuettel.com/code/rcpp/Rcpp-FAQ.pdf)

## Example

``` r
library(gpss)
#> 
#> Attaching package: 'gpss'
#> The following object is masked from 'package:stats':
#> 
#>     kernel
data(lalonde)
cat_vars=c("race_ethnicity", "married")
all_vars= c("age","educ","re74","re75","married", "race_ethnicity")

X = lalonde[,all_vars]


Y = lalonde$re78
D = lalonde$nsw

X_train = X[D==0,]
Y_train = Y[D==0]

X_test <- X[D == 1,]
Y_test <- Y[D == 1]

gp_train.out <- gp_train(X = X_train, Y = Y_train, optimize=TRUE,
                         mixed_data = TRUE, 
                         cat_columns = cat_vars)

gp_predict.out <- gp_predict(gp_train.out, X_test)
```
