\name{crq}
\alias{crq}
\alias{print.crq}
\title{Function to fit censored quantile regression model} 
\description{Fits a conditional quantile regression model for censored data.}
\usage{
crq(formula, data, subset, weights, na.action, method = "grid", contrasts = NULL, ...)
}
\arguments{
  \item{formula}{A formula object, with the response on the left of the `~'
          operator, and the terms on the right.  The response must be a
          survival object as returned by the \code{Surv} function.}
  \item{data}{A data.frame in which to interpret the variables named in the
          `formula',  in the `subset', and the `weights' argument.}
  \item{subset}{an optional vector specifying a subset of observations to be
          used in the fitting process.}
  \item{weights}{vector of observation weights; if supplied, the algorithm
          fits to minimize the sum of the weights multiplied into the
          absolute residuals. The length of weights vector must be the same as
          the number of observations.  The weights must be nonnegative
          and it is strongly recommended that they be strictly
          positive, since zero weights are ambiguous.} 
  \item{na.action}{a function to filter missing data.  This is applied to the
          model.frame after any subset argument has been used.  The
          default (with 'na.fail') is to create an error if any missing
          values are   found.  A possible alternative is 'na.omit',
          which  deletes observations that contain one or more missing
          values. } 
  \item{method}{The algorithmic method used for fitting.  There are currently
	two options.  The default method is "grid" which fits the crq process
	on a specified grid of tau values in (0,1), with default spacing
	\eqn{=\min (0.01,1/(2n^{.7}))}, where \eqn{n} is the sample size.
	The other method is "pivot" which fits the process on a grid
	determined by the data; this method is analogous to fitting the
	whole rq process with \code{rq}.}
  \item{contrasts}{a list giving contrasts for some or all of the factors 
          default = 'NULL' appearing in the model formula.  The
          elements of the list should have the same name as the
          variable  and should be either a contrast matrix
          (specifically, any full-rank  matrix with as many rows as
          there are levels in the factor),  or else a function to
          compute such a matrix given the number of levels.} 
  \item{...}{additional arguments for the fitting routine  (see
          \code{crq.fit}. } 
}
\details{ The algorithm uses recursively reweighted quantile regression to estimate 
       linear conditional quantile functions  for censored data as in Portnoy (2003).
       If method is "pivot", all break points are computed.  If method is "grid",
       then estimates are computed only on a fixed grid, see \code{crq.fit} for
	details on how to specify this grid.  Control parameters needed by
	\code{crq.fit} can be passed in \code{crq} via the \code{...} argument.
       For smaller data sets (n<400 or 500), using "pivot" is generally better
       and sufficiently quick. For n>400, it is generally adequate to 
       use "grid"  and let step be somewhat less that 1/sqrt(n).}
\value{see \code{\link{crq.object}}}
}
\references{Portnoy, Stephen (2003) Censored Regression Quantiles, JASA, 98, 1001-1010.} 

\author{Stephen Portnoy & Tereza Neocleous}

\seealso{\code{\link{crq.fit}},\code{\link{crq.object}}}

\examples{

data(uis)

#estimate the crq model using log(TIME) AFT specification
fit <- crq(Surv(log(TIME), CENSOR) ~  ND1 + ND2 + IV3 +
               TREAT + FRAC + RACE + AGE * SITE, data = uis)
Sfit <- summary(fit,1:19/20)
PHit <- coxph(Surv(TIME, CENSOR) ~  ND1 + ND2 + IV3 +
               TREAT + FRAC + RACE + AGE * SITE, data = uis)
plot(Sfit, CoxPHit = PHit)
}
\keyword{survival}
