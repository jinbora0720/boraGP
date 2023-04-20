\name{Cov_matern}
\alias{Cov_matern}
\title{Matern covariance matrix}

\description{
Compute a Matern covariance matrix
}

\usage{
Cov_matern(dist, sigmasq, phi, nu)
}

\arguments{
  \item{dist}{an \eqn{n \times n}{n x n} matrix of Euclidean distance.}
  
  \item{sigmasq}{partial sill parameter}
  
  \item{phi}{spatial decay parameter}
    
  \item{nu}{smoothness parameter}
}

\value{
  The return object is an \eqn{n \times n}{n x n} covariance matrix 
  corresponding to \code{dist}. Each element is computed by an isotropic 
  Matern covariance function: 
  \deqn{C(d) = \frac{\sigma^2(\phi d)^{\nu}}{\Gamma(\nu)2^{\nu-1}}
  K_{\nu}(\phi d)}
  where \eqn{K_{\nu}} is the modified Bessel function of the second kind 
  of order \eqn{\nu}.
}


\references{
  Jin, B., Herring, A.H., Dunson, D. (2022)
  Spatial predictions on physically constrained domains: 
  Applications to Arctic sea salinity data. 
  \href{https://arxiv.org/abs/2210.03913}{arXiv:2210.03913}.
}

\author{
  Bora Jin \email{bora.jin@duke.edu}
}

\examples{
  # create coordinates
  coords <- expand.grid(seq(0, 3, length = 10), seq(0, 3, length = 10))
  
  # compute Euclidean distance and convert it to a matrix
  dist <- as.matrix(dist(coords))
  
  # compute matern covariance matrix
  cov <- Cov_matern(dist = dist, sigmasq = 1, phi = 0.5, nu = 1.5)
}

\keyword{covariance}
