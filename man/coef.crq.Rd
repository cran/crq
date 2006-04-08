\name{coef.crq}
\alias{coef.crq}
\title{Extract coefficients from a crq object }
\description{This function extracts linearly interpolated coefficients
at specified quantiles from an object of class crq.}

\usage{coef.crq(object, taus = 1:4/5, ...) }
\arguments{
  \item{object}{An object of class "crq" produced by \code{crq}}
  \item{taus}{A vector of quantiles at which to evaluate the coefficients} 
  \item{...}{Other arguments passed for printing, etc.}
}
\details{Extracts coefficients by linear interpolation.}
\value{A matrix with the estimated coefficients at each of the chosen quantiles. 
The last row is the conditional quantile function evaluated at the mean design
point, i.e. Qbar = xbar ' betahat (tau).}

\author{Stephen Portnoy & Tereza Neocleous}

\seealso{\code{\link{crq}} }

\examples{ }

\keyword{survival}
