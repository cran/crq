\name{plot.summary.crqs}
\alias{plot.summary.crqs}
\title{Function to plot the effects from a censored regression quantile model} 
\description{This function plots the quantile treatment effects 
from a crq fit. It also (optionally) plots the corresponding Cox model 
estimates and a bootstrap confidence band.} 
\usage{
plot.summary.crqs(x,  nrow = 3, ncol = 3,  CoxPHit = NULL, ...)
}
\arguments{
  \item{x}{An object of class summary.crqs produced by \code{summary.crq}.} 
  \item{CoxPHit}{An object of class coxph produced by \code{coxph}, by
	default NULL indicating that no Cox effects are plotted.} 
  \item{nrow}{number of rows in plot array.}
  \item{ncol}{number of columns in plot array.}
  \item{...}{other arguments to the plotting function}
}

\details{The betahat (tau) effects estimated by the function \code{crq}
	are plotted, together with a confidence band produced by the
	function \code{summary.crq}.  Optionally, one can also overplot
	the corresponding estimates of the same covariate effects from
	a fitted Cox proportional hazard model.  See \code{CoxQTE}
	for details on how the Cox estimates are translated into 
	quantile treatment effects comparable to the crq estimates.
	There is no formal checking that the models implicit in the
	summary.crqs and coxph objects are the same.}


\seealso{\code{\link{crq}},\code{\link{summary.crq}},\code{\link{QTECox}}}

\keyword{survival}
