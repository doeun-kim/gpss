% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gp_functions.R
\name{gp_train}
\alias{gp_train}
\title{gp_train}
\usage{
gp_train(
  X,
  Y,
  b = NULL,
  s2 = 0.3,
  optimize = FALSE,
  scale = TRUE,
  mixed_data = FALSE,
  cat_columns = NULL,
  Xtest = NULL
)
}
\arguments{
\item{X}{a set of covariate data frame or matrix}

\item{Y}{Y vector (outcome variable)}

\item{b}{bandwidth (default = NULL)}

\item{s2}{noise or a fraction of Y not explained by X (default = 0.3)}

\item{optimize}{a logical value to indicate whether an automatic optimized value of S2 should be used. If FALSE, users must define s2. (default = FALSE)}

\item{scale}{a logical value to indicate whether covariates should be scaled. (dafault = TRUE)}

\item{mixed_data}{a logical value to indicate whether the covariates contain a categorical/binary variable (default = FALSE)}

\item{cat_columns}{a character or a numerical vector indicating categorical variables (default = NULL)}

\item{Xtest}{a data frame or a matrix of testing covariates. This is necessary when a non-overlapping categorical value exists between training and testing data sets. (default = NULL)}
}
\value{
\item{post_mean_scaled}{posterior distribution of Y in a scaled form}
\item{post_mean_orig}{posterior distribution of Y in an original scale}
\item{post_cov_scaled}{posterior covariance matrix in a scaled form}
\item{post_cov_orig}{posterior covariance matrix in an original scale}
\item{K}{a kernel matrix of X}
\item{prior_mean_scaled}{prior distribution of mean in a scaled form}
\item{X.orig}{the original matrix or data set of X}
\item{X.init}{the original matrix or data set of X with categorical variables in an expanded form}
\item{X.init.mean}{the initial mean values of X}
\item{X.init.sd}{the initial standard deviation values of X}
\item{Y.init.mean}{the initial mean value of Y}
\item{Y.init.sd}{the initial standard deviation value of Y}
\item{K}{the kernel matrix of X}
\item{Y}{scaled Y}
\item{X}{scaled X}
\item{b}{bandwidth}
\item{s2}{sigma squared}
\item{alpha}{alpha value in Rasmussen and Williams (2006) p.19}
\item{L}{L value in Rasmussen and Williams (2006) p.19}
\item{mixed_data}{a logical value indicating whether X contains a categorical/binary variable}
\item{cat_columns}{a character or a numerical vector indicating the location of categorical/binary variables in X}
\item{cat_num}{a numerical vector indicating the location of categorical/binary variables in an expanded version of X}
\item{Xcolnames}{column names of X}
}
\description{
to train GP model with training data set
}
\examples{
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

gp_train.out <- gp_train(X = X_train, Y = Y_train, optimize=TRUE, mixed_data = TRUE, cat_columns = cat_vars)
gp_predict.out <- gp_predict(gp_train.out, X_test)
}
