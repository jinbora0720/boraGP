\name{barrier_neighbor}
\alias{barrier_neighbor}
\title{BORA-GP neighbor search}

\description{
Find physically sensible nearest neighbors in a given domain
}

\usage{
barrier_neighbor(coords, coords.0 = NULL, ord = NULL,  
                 n.neighbors = 15, barrier, 
                 cores = 2, 
                 verbose = TRUE, 
                 debug = list(barrier_n.indx = NULL, 
                              barrier_dist = NULL, 
                              barrier_nn.indx.0_list = NULL, 
                              barrier_dist0 = NULL, 
                              ref_window = NULL, 
                              nonref_window = NULL,
                              ref_fill = TRUE, 
                              nonref_fill = TRUE))
}

\arguments{
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of reference coordinates 
  in \eqn{R^2}{R^2} (e.g., easting and northing). 
  There can be no duplicate locations.}
  
  \item{coords.0}{a matrix of non-reference coordinates 
  in \eqn{R^2}{R^2} (e.g., easting and northing). 
  There can be no duplicate locations.}
  
  \item{ord}{an index vector of length \eqn{n}{n} used for the neighbor search. 
  The default used when \code{ord} is not specified is x-axis ordering, i.e.,
  \code{order(coords[,1])}.}
    
  \item{n.neighbors}{number of neighbors used in BORA-GP.}
  
  \item{barrier}{an sf object of a given domain}
  
  \item{cores}{number of cores to use in parallel computing.}
  
  \item{verbose}{boolean to print algorithm progress.}
  
  \item{debug}{list for debugging. 
  Second-order neighbor search is activated for reference locations and non-reference 
  locations by \code{ref_fill = TRUE} and \code{nonref_fill = TRUE}, respectively. 
  With zero first-order neighbors, an imaginary grid is applied and users can 
  increase or decrease the number of locations to examine if crossing at each grid point using
  \code{ref_window} and \code{nonref_window}. 
  They are integer vectors of length \eqn{n}{n} and \eqn{nrow(coords.0)}{nrow(coords.0)},
  respectively. 
  If unspecified, the window size is fixed at \code{n.neighbors}.
  }
}

\value{
  The return object is a list comprising:
  \item{barrier_n.indx}{a list of BORA-GP neighbors for \code{coords}}
  \item{barrier_dist}{a list of distances to BORA-GP neighbors for \code{coords}}
  \item{barrier_nn.indx.0_list}{a list of BORA-GP neighbors for \code{coords.0}. 
  Returned if \code{coords.0} is specified.}
  \item{barrier_dist0}{a list of distances to BORA-GP neighbors for \code{coords.0}. 
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
