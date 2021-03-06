% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crrHARD.R
\name{crrHARD}
\alias{crrHARD}
\title{Hard Threshold Regression for Competing Risks Regression}
\usage{
crrHARD(ftime, fstatus, X, failcode = 1, cencode = 0, lambda = 0,
  nlambda = 25)
}
\arguments{
\item{ftime}{A vector of event/censoring times.}

\item{fstatus}{A vector with unique code for each event type and a separate code for censored observations.}

\item{X}{A matrix of fixed covariates (nobs x ncovs)}

\item{failcode}{Integer: code of \code{fstatus} that event type of interest (default is 1)}

\item{cencode}{Integer: code of \code{fstatus} that denotes censored observations (default is 0)}

\item{nlambda}{Numeric: number of \code{lambda} values if performing grid search  (default is 25)}
}
\value{
Returns a list of class \code{crrBAR}.
}
\description{
Fits hard thresholding penalty for competing risks regression. Estimates are calculated
via \code{crr} from the \strong{cmprsk} package.
}
\details{
The \code{crrHARD} performs hard thresholding on the MPLE using the rule \eqn{\hat{\beta_j} = \hat{beta}_j I(|\beta_j| > \lambda)}
}
\examples{
set.seed(10)
ftime <- rexp(200)
fstatus <- sample(0:2, 200, replace = TRUE)
cov <- matrix(runif(1000), nrow = 200)
dimnames(cov)[[2]] <- c('x1','x2','x3','x4','x5')
fit <- crrHARD(ftime, fstatus, cov, lambda = 0.05)
fit$coef
}
\references{
Fine J. and Gray R. (1999) A proportional hazards model for the subdistribution of a competing risk.  \emph{JASA} 94:496-509.
}
