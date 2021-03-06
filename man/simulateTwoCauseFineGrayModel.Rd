% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulateTwoCauseFineGrayModel.R
\name{simulateTwoCauseFineGrayModel}
\alias{simulateTwoCauseFineGrayModel}
\title{Simulate data from the Fine-Gray Model}
\usage{
simulateTwoCauseFineGrayModel(nobs, beta1, beta2, X = NULL, u.min = 0,
  u.max, p = 0.5, returnX = FALSE)
}
\arguments{
\item{nobs}{Integer: Number of observations in simulated dataset.}

\item{beta1}{A vector of effect sizes for cause 1 of length ncovs}

\item{beta2}{A vector of effect sizes for cause 2 of length ncovs}

\item{X}{A matrix of fixed covariates (nobs x ncovs). If \code{X} is NULL (default) then \code{X} will be simulated from MVN(O, I) with n = nobs and p = length(beta1).}

\item{u.min}{Numeric: controls lower bound of censoring distribution where C ~ U(u.min, u.max)}

\item{u.max}{Numeric: controls upper bound of censoring distribution where C ~ U(u.min, u.max)}

\item{p}{Numeric: value between 0 and 1 which controls the mixture probability.}

\item{returnX}{Logical: Whether to return \code{X} or not. Default is TRUE. Recommended if \code{X} is NULL.}
}
\value{
Returns a list with \code{ftime}, \code{fstatus}, \code{X}.
}
\description{
Simulate data from the model proposed in Fine and Gray (1999) for two causes. Cause 1 is assumed
to be of primary importance.
}
\details{
The function simulates data according to the setup by Fine and Gray (1999). See their paper for more information.
}
\examples{
set.seed(10)
nobs <- 500
beta1 <- c(0.5, 0.4, 0, 0, 0.35, 0, 0, 0.7)
beta2 <- -beta1
X <- matrix(rnorm(nobs * 8), nrow = nobs)
dat <- simulateTwoCauseFineGrayModel(nobs, beta1, beta2, X, u.min = 0, u.max = 1, p = 0.5)
}
\references{
Fine J. and Gray R. (1999) A proportional hazards model for the subdistribution of a competing risk.  \emph{JASA} 94:496-509.
}
