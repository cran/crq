\name{summary.crq}
\alias{summary.crq}
\alias{print.summary.crq}
\alias{print.summary.crqs}
\title{Summary function for censored quantile regression estimation}

\description{Confidence intervals and standard error computations for
censored quantile regression.}
\usage{
summary.crq(object, taus = 1:4/5, alpha = .1, se = "boot", ...)
}
\arguments{
  \item{object}{An object of class \code{crq}.}
  \item{taus}{A vector of quantiles at which to evaluate estimates, their
	confidence intervals and standard errors.}
  \item{alpha} {Uncoverage probability for the reported confidence intervals}
  \item{se}{A method for estimating standard errors.}
  \item{...}{Other arguments to be passed to \code{boot.crq}, etc.}
}
\details{Only one method is (so far) available for confidence interval
and standard error computations; it uses the bootstrap to compute rescaled
interquartile range estimates.  Two flavors of the bootstrap are available:
by default one obtains the conventional xy-pair bootstrap, if bmethod = "Bose"
a weighted resampling method with standard exponential weights is used.
Confidence intervals are computed by rescaling the interquartile range,
of the bootstrap realizations as are the standard errors as suggested in Portnoy (2003).}
\value{A list of class summary.crqs each member consists of a tau component
and a coefficient component.  The latter component is a matrix containing
the point, interval and standard error estimates.}

\author{Stephen Portnoy & Tereza Neocleous}
\seealso{\code{boot.crq}, \code{crq}}

\examples{ }
\keyword{survival }
