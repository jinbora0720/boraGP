\name{rboraGP}
\alias{rboraGP}
\title{Sample data from BORA-GP}

\description{
Generate random data from prior sampling of univariate Bayesian
regression models based on BORA-GP with Gaussian errors
}

\usage{
rboraGP(coords, neighbor.info, X = NULL,
        beta, sig_sq, phi, nu = NULL, tau_sq,
        coords.0 = NULL, nn.indx.0, X.0 = NULL,
        base_cov = "exponential")
}

\arguments{
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of reference coordinates 
  in \eqn{R^2}{R^2} (e.g., easting and northing). 
  There can be no duplicate locations.}
  
  \item{neighbor.info}{a list of \code{ord} and \code{n.indx}. 
  They are outputs of \code{\link{barrier_neighbor}}.}
  
  \item{X}{an \eqn{n \times p}{n x p} matrix of covariates.}
  
  \item{beta}{a vector of regression coefficients of length \code{p}.}
  
  \item{sig_sq}{a partial sill.}
  
  \item{phi}{a spatial decay.}
  
  \item{nu}{a smoothness parameter, applicable only for \code{base_cov = "matern"}.}
  
  \item{tau_sq}{a nugget.}
  
  \item{coords.0}{a matrix of non-reference coordinates 
  in \eqn{R^2}{R^2} (e.g., easting and northing). 
  There can be no duplicate locations.}
  
  \item{nn.indx.0}{a neighbor info matrix for coords.0 of size 
  \eqn{nrow(coords.0) \times m}{nrow(coords.0) x m} 
  where \eqn{m}{m} is the number of neighbors.}
  
  \item{X.0}{a matrix of covariates for \code{coords.0}.}

  \item{base_cov}{a base covariance function. Either 
  \code{"exponential"} or \code{"matern"} are supported.}  
}

\value{

  The return object is a list comprising:
  \item{y}{sampled outcome \code{y} for \code{coords}.}
  \item{w}{sampled latent process \code{w} for \code{coords}.}
  \item{y.0}{sampled outcome \code{y} for \code{coords.0}. 
  Returned if \code{coords.0} is specified.}
  \item{w.0}{sampled latent process \code{w} for \code{coords.0}. 
  Returned if \code{coords.0} is specified.}

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
}

\keyword{model}
