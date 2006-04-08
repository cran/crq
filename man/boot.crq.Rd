\name{boot.crq}
\alias{boot.crq}
\title{Bootstrap for censored regresion quantile models} 
\description{This function returns an array of bootstrap estimates for censored 
regression quantile models. }
\usage{ boot.crq(x, y, c, taus, R = 100, mboot, bmethod = "xy-pair", ...) }
\arguments{
  \item{x}{The design matrix}
  \item{y}{The response vector}
  \item{c}{The censoring indicator}
  \item{taus}{The vector of taus at which to evaluate.}
  \item{R}{Number of bootstrap replications, defaults to 100.}
  \item{mboot}{Subsample size, defaults to the original sample size.}
  \item{bmethod}{method for the bootstrap, either "xy-pair" or "Bose".}
  \item{...}{Other parameters to be passed to \code{crq.fit}.}
}
\details{Resampling function for crq.  If method = "xy-pair", the default,
	then R samples of size mboot are taken and crq estimates are computed. 
	If method = "Bose" the generalized resampling strategy of drawing
	standard exponential weights is used.}

\value{ Returns an list of three components, 
	the of dimension p+1 by length(taus) by R, where
	p is the row dimension of the design matrix, i.e. the dimension
	of the model parameter.}

\author{Stephen Portnoy & Tereza Neocleous}

\seealso{\code{\link{crq}} }

\keyword{survival}
