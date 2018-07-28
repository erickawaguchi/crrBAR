#' Hard Threshold Regression for Competing Risks Regression
#'
#' @description Fits hard thresholding penalty for competing risks regression. Estimates are calculated
#' via \code{crr} from the \strong{cmprsk} package.
#'
#' @param ftime A vector of event/censoring times.
#' @param fstatus A vector with unique code for each event type and a separate code for censored observations.
#' @param X A matrix of fixed covariates (nobs x ncovs)
#' @param failcode Integer: code of \code{fstatus} that event type of interest (default is 1)
#' @param cencode Integer: code of \code{fstatus} that denotes censored observations (default is 0)
#' @param nlambda Numeric: number of \code{lambda} values if performing grid search  (default is 25)
#'
#' @details The \code{crrHARD} performs hard thresholding on the MPLE using the rule \eqn{\hat{\beta_j} = \hat{beta}_j I(|\beta_j| > \lambda)}
#'
#' @return Returns a list of class \code{crrBAR}.
#'
#' @import survival cmprsk
#' @export
#' @useDynLib crrBAR standardize
#' @examples
#' set.seed(10)
#' ftime <- rexp(200)
#' fstatus <- sample(0:2, 200, replace = TRUE)
#' cov <- matrix(runif(1000), nrow = 200)
#' dimnames(cov)[[2]] <- c('x1','x2','x3','x4','x5')
#' fit <- crrHARD(ftime, fstatus, cov, lambda = 0.05)
#' fit$coef
#' @references
#' Fine J. and Gray R. (1999) A proportional hazards model for the subdistribution of a competing risk.  \emph{JASA} 94:496-509.

crrHARD <- function(ftime, fstatus, X, failcode = 1, cencode = 0,
                   lambda = 0,
                   nlambda = 25){

  # Function to calculate generalized inverse
  ginv = function(X, tol = sqrt(.Machine$double.eps)){
    s = svd(X)
    nz = s$d > tol * s$d[1]
    if(any(nz)) s$v[,nz] %*% (t(s$u[, nz]) / s$d[nz])
    else X*0
  }

  #- Calculate MPLE via crr
  # Standardize design matrix here
  std    <- .Call("standardize", X, PACKAGE = "crrBAR")
  XX     <- std[[1]]
  center <- std[[2]]
  scale  <- std[[3]]
  nz <- which(scale > 1e-6)
  if (length(nz) != ncol(XX)) XX <- XX[ , nz, drop = FALSE]

  fit.crr <- crr(ftime, fstatus, cov1 = XX, failcode = failcode,
                 cencode = cencode)

  bhat <- fit.crr$coef
  # If lambda is MISSING, create lambda path
  if(missing(lambda)) {
    lambda <- seq(0, max(abs(bhat)), length = nlambda)
  } else {
    lambda <- lambda
  }
  if(min(lambda) < 0) stop("lambda must be a non-negative number.")

  #Results to store:
  coefMatrix           <- matrix(NA, nrow = ncol(X), ncol = length(lambda))
  colnames(coefMatrix) <- round(lambda, 3)

  logLik  <- numeric(length(lambda))

  #Enter BAR Fit here (for different lambdas)
  for(l in 1:length(lambda)) {
    lam <- lambda[l]
    betaHT <- ifelse(abs(bhat) > lam, bhat, 0) / scale
    coefMatrix[, l] <- betaHT
    logLik[l] <- getLogLikelihood(ftime, fstatus, X, failcode = 1,
                                  cencode = 0, beta = betaHT)
  }

  ## Output
  val <- structure(list(coef = coefMatrix,
                        logLik = logLik,
                        call = sys.call()),
                   class = "crrBAR")

  val
}
