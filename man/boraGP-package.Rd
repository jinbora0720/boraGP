\name{boraGP-package}
\alias{boraGP-package}
\alias{boraGP}
\docType{package}
\title{
\packageTitle{boraGP}
}
\description{
Fits univariate Bayesian spatial regressions for large data using Barrier Overlap-Removal Acyclic Directed Graph Gaussian Process (BORA-GP). 
}
\details{
The function \code{\link{barrier_neighbor}} finds physically sensible nearest 
neighbors in a given domain. 

The function \code{\link{spNNGP}} conducts posterior sampling using
Markov chain Monte Carlo for univariate Bayesian regression models. For BORA-GP, 
one can specify neighbors as the outcome of \code{\link{barrier_neighbor}}. If 
neighbors are not specified or barriers do not exist, it comes down to a usual 
Nearest Neighbor GP (NNGP; Datta et al. (2016)).

The function \code{\link{predict}} makes predictions on new locations based on 
estimation results by \code{\link{spNNGP}}. 

The function \code{\link{rboraGP}} generates random data from prior sampling of
univariate Bayesian regression models based on BORA-GP with Gaussian errors.
}
\author{
Bora Jin \href{mailto:bora.jin@duke.edu}{bora.jin@duke.edu}
}
\references{
Jin, B., Herring, A.H., Dunson, D. (2022). Scalable Gaussian Processes on Physically Constrained Domains. \href{}{arXiv}

The code heavily depends on Finley, Datta, Banerjee (2022) \href{https://doi.org/10.18637/jss.v103.i05}{doi:10.18637/jss.v103.i05}, Finley, Datta, Cook, Morton, Andersen, and Banerjee (2019) \href{https://doi.org/10.1080/10618600.2018.1537924}{doi:10.1080/10618600.2018.1537924} and Datta, Banerjee, Finley, and Gelfand (2016) \href{https://doi.org/10.1080/01621459.2015.1044091}{doi:10.1080/01621459.2015.1044091}.
}
\keyword{ package }
\seealso{
\code{\link{rboraGP}}, \code{\link{barrier_neighbor}, \code{\link{spNNGP}}, \code{\link{predict.NNGP}}}
}
