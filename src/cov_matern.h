#include <RcppArmadillo.h>
#include <R.h>

arma::mat Cov_matern(const arma::mat& dist,
                     const double& sigmasq,
                     const double& phi, const double& nu);
