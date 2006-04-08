\name{crq.fit}
\alias{crq.fit}
\title{Fitting of censored quantile regression models} 

\description{Function controlling the fitting of censored quantile regression models} 
\usage{ crq.fit(x, y, cen, weights, method, mw = 20 , ginit = 0.001, gstep = NULL) }
\arguments{
  \item{x}{the design  matrix}
  \item{y}{the response vector typically in log(event time) (AFT) form}
  \item{cen}{ censoring indicator:  usually  0 for censored, 
	1 for uncensored. No provision for interval censoring (yet).}
  \item{weights}{weights to be used for fitting see \code{crq}.}
  \item{method}{method to be used: either "grid" or "pivot".}
  \item{mw}{Maximal number of allowed calls to weighted rq 
	(to handle possible degeneracy).  Only relevant for method = "pivot".}
   \item{gstep}{spacing  of points of evaluation only relevant for method = "grid",
        defaults to \code{min(0.01,1/(2*length(y)^{.7}))}.} 
   \item{ginit}{initial point of evaluation only relevant for method = "grid",
        defaults to 0.001.} 
   }
\details{See {\code{crq}}} 
\value{
  \item{sol}{A matrix with (p+2) rows:  the first row is the taus (in (0,1))
	corresponding to the fitted model parameters.  The next p rows are
	the fitted parameters betahat (tau), and the final row contains 
	entries predicting the conditional quantiles of the response at
	the mean of the design:  xbar ' betahat (tau).}
  \item{Isplit}{Indices of split censored observations.} 
  \item{tausplit}{Corresponding tau values.}
  \item{status}{Status of censored observations: 
            0: not censored 
            1: split (crossed) censored 
            2: deleted censored as below tau = 0 solution 
            3: above last (maximal tau) solution}
}

\author{Stephen Portnoy & Tereza Neocleous}
\seealso{\code{\link{crq}}}

\keyword{survival}
