#include <RcppArmadillo.h>
#include <R.h>
//#include "cov_matern.h"

//[[Rcpp::export]]
arma::mat Cov_matern(const arma::mat& dist, 
                     const double& sigmasq, const double& phi, const double& nu){
  
  int n = dist.n_rows;
  arma::mat out = arma::zeros<arma::mat>(n, n);
  double pow2_nu1_gammanu = pow(2.0, 1.0-nu) / R::gammafn(nu);
  for (int i=0; i < n; i++) {
    for (int j=i+1; j < n; j++) {
      double hphi = dist(i,j)*phi;
      if (hphi > 0.0) {
        out(i, j) = sigmasq * pow(hphi, nu) * pow2_nu1_gammanu *
          R::bessel_k(hphi, nu, 1.0);
      } else {
        out(i, j) = sigmasq;
      }
    }
  }
  out += out.t();
  out.diag().fill(sigmasq);
  
  return out;
}



