#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// 1. Gaussian kernel
// [[Rcpp::export]]
arma::mat kernel_gaussian(
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
arma::mat kernel_symmetric_gaussian(
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

// 2. Gaussian + linear
// [[Rcpp::export]]
arma::mat kernel_gaussian_linear(
  arma::mat x1,
  arma::mat x2,
  double b
) {
  size_t n1 = x1.n_rows;
  size_t n2 = x2.n_rows;
  arma::mat K(n1,n2);
  double sqdist, gaussian_part, linear_part;
  size_t i,j;
  
  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++) {
      // Gaussian part
      sqdist = arma::sum(square(x1.row(i)-x2.row(j)));
      gaussian_part = std::exp(-sqdist/b);
      
      // Linear part
      linear_part = arma::dot(x1.row(i), x2.row(j));
      
      K(i,j) = gaussian_part + linear_part;
    }
  }
  return(K);
}

// [[Rcpp::export]]
arma::mat kernel_symmetric_gaussian_linear(
  arma::mat x,
  double b
) {
  size_t n = x.n_rows;
  arma::mat K(n, n);
  double sqdist, gaussian_part, linear_part;
  
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i; j < n; j++) {
      if (i == j) {
        // Diagonal: gaussian = 1, plus linear self-product
        K(i,j) = 1.0 + arma::dot(x.row(i), x.row(i));
      } else {
        // Off-diagonal
        sqdist = arma::sum(arma::square(x.row(i) - x.row(j)));
        gaussian_part = std::exp(-sqdist / b);
        linear_part = arma::dot(x.row(i), x.row(j));
        K(i,j) = gaussian_part + linear_part;
        K(j,i) = K(i,j); // Symmetric property
      }
    }
  }
  return K;
}

// 3. Gaussian + quadratic
// [[Rcpp::export]]
arma::mat kernel_gaussian_quadratic(
  arma::mat x1,
  arma::mat x2,
  double b
) {
  size_t n1 = x1.n_rows;
  size_t n2 = x2.n_rows;
  arma::mat K(n1,n2);
  double sqdist, gaussian_part, quadratic_part, dot_product;
  size_t i,j;
  
  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++) {
      // Gaussian part
      sqdist = arma::sum(square(x1.row(i)-x2.row(j)));
      gaussian_part = std::exp(-sqdist/b);
      
      // Quadratic part
      dot_product = arma::dot(x1.row(i), x2.row(j));
      quadratic_part = std::pow(dot_product + 1.0, 2.0);
      
      K(i,j) = gaussian_part + quadratic_part;
    }
  }
  return(K);
}

// [[Rcpp::export]]
arma::mat kernel_symmetric_gaussian_quadratic(
  arma::mat x,
  double b
) {
  size_t n = x.n_rows;
  arma::mat K(n, n);
  double sqdist, gaussian_part, quadratic_part, dot_product;
  
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i; j < n; j++) {
      dot_product = arma::dot(x.row(i), x.row(j));
      
      if (i == j) {
        // Diagonal: gaussian = 1, plus quadratic
        quadratic_part = std::pow(dot_product + 1.0, 2.0);
        K(i,j) = 1.0 + quadratic_part;
      } else {
        // Off-diagonal
        sqdist = arma::sum(arma::square(x.row(i) - x.row(j)));
        gaussian_part = std::exp(-sqdist / b);
        quadratic_part = std::pow(dot_product + 1.0, 2.0);
        K(i,j) = gaussian_part + quadratic_part;
        K(j,i) = K(i,j); // Symmetric property
      }
    }
  }
  return K;
}

// 4. Gaussian + periodic + linear
// [[Rcpp::export]]
arma::mat kernel_gaussian_periodic_linear(
  arma::mat x1,
  arma::mat x2,
  double b,
  double period
) {
  size_t n1 = x1.n_rows;
  size_t n2 = x2.n_rows;
  size_t ncol = x1.n_cols;
  arma::mat K(n1,n2);
  double sqdist, gaussian_part, periodic_sum, linear_part, diff;
  size_t i,j,d;
  
  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++) {
      // Gaussian part
      sqdist = arma::sum(square(x1.row(i)-x2.row(j)));
      gaussian_part = std::exp(-sqdist/b);

      // Periodic part with b/2 as lengthscale
      double periodic_sum = 0.0;
      for (d=0; d<ncol; d++) {
        diff = x1(i,d) - x2(j,d);
        periodic_sum += 2.0 * std::pow(std::sin(M_PI * std::abs(diff) / period), 2.0);
      }
      double periodic_part = std::exp(-periodic_sum / (b/2.0));
      
      // Linear part
      linear_part = arma::dot(x1.row(i), x2.row(j));
      
      K(i,j) = gaussian_part + periodic_part + linear_part;
    }
  }
  return(K);
}

// [[Rcpp::export]]
arma::mat kernel_symmetric_gaussian_periodic_linear(
  arma::mat x,
  double b,
  double period
) {
  size_t n = x.n_rows;
  size_t ncol = x.n_cols;
  arma::mat K(n, n);
  double sqdist, gaussian_part, periodic_sum, linear_part, diff;
  
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i; j < n; j++) {
      if (i == j) {
        // Diagonal: gaussian = 1, periodic = 1 (sin(0)=0, so exp(0)=1), plus linear self-product
        linear_part = arma::dot(x.row(i), x.row(i));
        K(i,j) = 1.0 + 1.0 + linear_part;
      } else {
        // Off-diagonal
        // Gaussian part
        sqdist = arma::sum(arma::square(x.row(i) - x.row(j)));
        gaussian_part = std::exp(-sqdist / b);

        // Periodic part with b/2 as lengthscale
        double periodic_sum = 0.0;
        for (size_t d = 0; d < ncol; d++) {
          diff = x(i,d) - x(j,d);
          periodic_sum += 2.0 * std::pow(std::sin(M_PI * std::abs(diff) / period), 2.0);
        }
        double periodic_part = std::exp(-periodic_sum / (b/2.0));
        
        // Linear part
        linear_part = arma::dot(x.row(i), x.row(j));
        
        K(i,j) = gaussian_part + periodic_part + linear_part;
        K(j,i) = K(i,j); // Symmetric property
      }
    }
  }
  return K;
}

// 5. Gaussian + periodic + quadratic
// [[Rcpp::export]]
arma::mat kernel_gaussian_periodic_quadratic(
  arma::mat x1,
  arma::mat x2,
  double b,
  double period
) {
  size_t n1 = x1.n_rows;
  size_t n2 = x2.n_rows;
  size_t ncol = x1.n_cols;
  arma::mat K(n1,n2);
  double sqdist, gaussian_part, periodic_sum, quadratic_part, dot_product, diff;
  size_t i,j,d;
  
  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++) {
      // Gaussian part
      sqdist = arma::sum(square(x1.row(i)-x2.row(j)));
      gaussian_part = std::exp(-sqdist/b);
      
      // Periodic part with b/2 as lengthscale
      double periodic_sum = 0.0;
      for (d=0; d<ncol; d++) {
        diff = x1(i,d) - x2(j,d);
        periodic_sum += 2.0 * std::pow(std::sin(M_PI * std::abs(diff) / period), 2.0);
      }
      double periodic_part = std::exp(-periodic_sum / (b/2.0));
      
      // Quadratic part
      dot_product = arma::dot(x1.row(i), x2.row(j));
      quadratic_part = std::pow(dot_product + 1.0, 2.0);
      
      K(i,j) = gaussian_part + periodic_part + quadratic_part;
    }
  }
  return(K);
}

// [[Rcpp::export]]
arma::mat kernel_symmetric_gaussian_periodic_quadratic(
  arma::mat x,
  double b,
  double period
) {
  size_t n = x.n_rows;
  size_t ncol = x.n_cols;
  arma::mat K(n, n);
  double sqdist, gaussian_part, periodic_sum, quadratic_part, dot_product, diff;
  
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i; j < n; j++) {
      dot_product = arma::dot(x.row(i), x.row(j));
      
      if (i == j) {
        // Diagonal: gaussian = 1, periodic = 1 (sin(0)=0, so exp(0)=1), plus linear self-product
        quadratic_part = std::pow(dot_product + 1.0, 2.0);
        K(i,j) = 1.0 + 1.0 + quadratic_part;
      } else {
        // Off-diagonal
        // Gaussian part
        sqdist = arma::sum(arma::square(x.row(i) - x.row(j)));
        gaussian_part = std::exp(-sqdist / b);
        
        // Periodic part with b/2 as lengthscale
        double periodic_sum = 0.0;
        for (size_t d = 0; d < ncol; d++) {
          diff = x(i,d) - x(j,d);
          periodic_sum += 2.0 * std::pow(std::sin(M_PI * std::abs(diff) / period), 2.0);
        }
        double periodic_part = std::exp(-periodic_sum / (b/2.0));
        
        // Quadratic part
        quadratic_part = std::pow(dot_product + 1.0, 2.0);
        
        K(i,j) = gaussian_part + periodic_part + quadratic_part;
        K(j,i) = K(i,j); // Symmetric property
      }
    }
  }
  return K;
}


// log likelihood
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



