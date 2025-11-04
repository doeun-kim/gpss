#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat kernel(
    arma::mat x1, // input matrix 1
    arma::mat x2, // input matrix 2
    double b // length-scale
) {
  size_t n1 = x1.n_rows;
  size_t n2 = x2.n_rows;
  arma::mat K(n1,n2);
  double sqdist;
  size_t i,j;

  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++) {
      sqdist = arma::sum(square(x1.row(i)-x2.row(j)));
      K(i,j) = std::exp(-sqdist/b);
    }
  }
  return(K);
}

// [[Rcpp::export]]
arma::mat kernel_symmetric(
    arma::mat x, // input matrix
    double b // length-scale
) {
  size_t n = x.n_rows;
  arma::mat K(n, n);
  double sqdist;

  K.diag().ones(); // Initialize the diagonal elements to 1

  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      sqdist = arma::sum(arma::square(x.row(i) - x.row(j)));
      K(i, j) = std::exp(-sqdist / b);
      K(j, i) = K(i, j); // Symmetric property
    }
  }

  return K;
}

// [[Rcpp::export]]
arma::mat kernel_linear_cpp(
    arma::mat x1, // input matrix 1
    arma::mat x2, // input matrix 2
    double sigma_f // magnitude
) {
  double sigma_f2 = sigma_f*sigma_f;
  arma::mat K = sigma_f2 * (x1 * x2.t());
  return(K);
}

// [[Rcpp::export]]
arma::mat kernel_se_cpp(
    arma::mat x1, // input matrix 1
    arma::mat x2, // input matrix 2
    double sigma_f, // magnitude
    double l // length-scale
) {
  size_t n1 = x1.n_rows;
  size_t n2 = x2.n_rows;
  arma::mat K(n1,n2);
  double l2 = l*l;
  double sigma_f2 = sigma_f*sigma_f;
  double sqdist;
  size_t i,j;

  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++) {
      sqdist = sum(square(x1.row(i)-x2.row(j)));
      K(i,j) = sigma_f2 * std::exp(-0.5*sqdist/l2);
    }
  }
  return(K);
}

// [[Rcpp::export]]
arma::mat kernel_periodic_cpp(
    const arma::mat& x1,
    const arma::mat& x2,
    double sigma_f, // magnitude
    double l, // length-scale
    double p // period
) {
  double l2 = l*l;
  double sigma_f2 = sigma_f*sigma_f;

  // Calculate the pairwise squared differences for each dimension
  arma::mat K(x1.n_rows, x2.n_rows);
  for (size_t d = 0; d < x1.n_cols; ++d) {
    arma::mat x1_d = repmat(x1.col(d), 1, x2.n_rows);
    arma::mat x2_d = repmat(x2.col(d).t(), x1.n_rows, 1);
    K += square(sin((x1_d - x2_d) * M_PI / p));
  }
  K = sigma_f2 * exp(-2 * K / l2);

  return K;
}


// [[Rcpp::export]]
double log_marginal_likelihood_cpp(
    const arma::mat& K,
    const arma::vec& y,
    double s2) {

  int n = K.n_rows;
  arma::mat L = arma::chol(K + s2 * arma::eye(n, n), "lower");
  arma::vec alpha = arma::solve(arma::trimatl(L), y, arma::solve_opts::fast + arma::solve_opts::no_approx);
  alpha = arma::solve(arma::trimatu(L.t()), alpha, arma::solve_opts::fast + arma::solve_opts::no_approx);

  double logDetK = arma::sum(arma::log(L.diag())); // Log determinant of K
  double logLik = -0.5 * arma::dot(y, alpha) - logDetK - (n / 2.0) * std::log(2 * M_PI);

  return(logLik);
}








