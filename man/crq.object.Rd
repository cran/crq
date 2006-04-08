\name{crq.object}
\alias{crq.object}
\title{Object of class `crq'}
\description{ This class of objects is returned by the `crq' class of
     functions to represent a fitted censored quantile regression model.}

\value{A model fit object containing inter alia:
\item{sol}{A matrix with p+2 rows:  the first row contains the values of
	tau, at which the crq process is evaluated (this is either the
	set of pivot values, or the set of grid points depending on the
	method argument of the call to \code{crq}).  
	The next p rows the estimated model parameters for the
	corresponding tau.  The last row contains the estimated conditional
	quantiles of the response evaluated at xbar.} 
\item{Isplit}{Indices of split censored observations.} 
\item{tausplit}{Corresponding tau values.}
\item{status}{Status of censored observations: 
            0: not censored, 
            1: split (crossed) censored,
            2: deleted censored as below tau = 0 solution, 
            3: above last (maximal tau) solution.} 
}

\seealso{\code{\link{crq}}}

\author{Stephen Portnoy & Tereza Neocleous}

\keyword{survival}
