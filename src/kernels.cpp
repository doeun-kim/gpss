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
  double g = 1e-8;
  bool success = false;
  arma::mat L;
  arma::mat I = arma::eye(n, n);
  
  while (!success && g < 1) {
    try {
      L = arma::chol(K + s2 * I + g * I, "lower");
      success = true;
    } catch(...) {
      g *= 10; // Increase g if chol fails
    }
  }
  
  if (!success) {
    stop("Cholesky decomposition failed. The kernel matrix might not be positive semi-definite.");
  }
  
  arma::vec alpha = arma::solve(arma::trimatl(L), y); // Forward substitution
  alpha = arma::solve(arma::trimatu(L.t()), alpha); // Backward substitution
  
  double logDetK = arma::as_scalar( arma::sum(arma::log(L.diag())) ); // Log determinant of K
  double logLik = arma::as_scalar( -0.5 * dot(y, alpha) - logDetK - n / 2.0 * std::log(2 * M_PI) );
  
  return(logLik);
}
