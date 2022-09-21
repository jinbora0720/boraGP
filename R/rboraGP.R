rboraGP <- function(coords, neighbor.info, X = NULL,
                    beta, sig_sq, phi, nu = NULL, tau_sq,
                    coords.0 = NULL, nn.indx.0, X.0 = NULL,
                    base_cov = "exponential") {

  Ctilde <- create_Ctilde(coords = coords,
                          neighbor.info = neighbor.info,
                          sig_sq = sig_sq, phi = phi, nu = nu,
                          base_cov = base_cov)
  
  n_S <- nrow(coords)
  w_S <- t(chol(Ctilde))%*%rnorm(n_S)
  if (!is.null(X)) {
    y_S <- X%*%beta + w_S + sqrt(tau_sq)*rnorm(n_S)
  } else {
    y_S <- w_S + sqrt(tau_sq)*rnorm(n_S)
  }
  
  if (!is.null(coords.0)) {
    n_U <- nrow(coords.0)
    w_U <- rep(0, n_U)
    if (base_cov == "exponential") {
      for (i in 1:n_U) {
        pidx <- nn.indx.0[i,]
        coords_tmp <- rbind(coords.0[i,], coords[pidx,])
        Cor <- exp(-phi*as.matrix(dist(coords_tmp)))
        Bu <- Cor[1,-1]%*%solve(Cor[-1,-1])
        Fu <- sig_sq*(Cor[1,1]-Bu%*%Cor[-1,1])
        w_U[i] <- Bu%*%w_S[pidx]+sqrt(Fu)*rnorm(1)
      }
    } else {
      for (i in 1:n_U) {
        pidx <- nn.indx.0[i,]
        coords_tmp <- rbind(coords.0[i,], coords[pidx,])
        Cor <- Cov_matern(as.matrix(dist(coords_tmp)), sigmasq = 1,
                          phi = phi, nu = nu)
        Bu <- Cor[1,-1]%*%solve(Cor[-1,-1])
        Fu <- sig_sq*(Cor[1,1]-Bu%*%Cor[-1,1])
        w_U[i] <- Bu%*%w_S[pidx]+sqrt(Fu)*rnorm(1)
      }
    }
    
    if (!is.null(X.0)) {
      y_U <- X.0%*%beta + w_U + sqrt(tau_sq)*rnorm(n_U)
    } else {
      y_U <- w_U + sqrt(tau_sq)*rnorm(n_U)
    }
    
    return(list(y=y_S, w=w_S, y.0=y_U, w.0=w_U))
  }
  
  return(list(y=y_S, w=w_S))
}