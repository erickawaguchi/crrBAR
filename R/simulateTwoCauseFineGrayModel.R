#' Simulate data from the Fine-Gray Model
#'
#' @description  Simulate data from the model proposed in Fine and Gray (1999) for two causes. Cause 1 is assumed
#' to be of primary importance.
#'
#' @param nobs Integer: Number of observations in simulated dataset.
#' @param beta1 A vector of effect sizes for cause 1 of length ncovs
#' @param beta2 A vector of effect sizes for cause 2 of length ncovs
#' @param X A matrix of fixed covariates (nobs x ncovs). If \code{X} is NULL (default) then \code{X} will be simulated from MVN(O, I) with n = nobs and p = length(beta1).
#' @param u.min Numeric: controls lower bound of censoring distribution where C ~ U(u.min, u.max)
#' @param u.max Numeric: controls upper bound of censoring distribution where C ~ U(u.min, u.max)
#' @param p Numeric: value between 0 and 1 which controls the mixture probability.
#' @param returnX Logical: Whether to return \code{X} or not. Default is TRUE. Recommended if \code{X} is NULL.
#' @details The function simulates data according to the setup by Fine and Gray (1999). See their paper for more information.
#' @return Returns a list with \code{ftime}, \code{fstatus}, \code{X}.
#' @examples
#' set.seed(10)
#' nobs <- 500
#' beta1 <- c(0.5, 0.4, 0, 0, 0.35, 0, 0, 0.7)
#' beta2 <- -beta1
#' X <- matrix(rnorm(nobs * 8), nrow = nobs)
#' dat <- simulateTwoCauseFineGrayModel(nobs, beta1, beta2, X, u.min = 0, u.max = 1, p = 0.5)
#' @references
#' Fine J. and Gray R. (1999) A proportional hazards model for the subdistribution of a competing risk.  \emph{JASA} 94:496-509.
#' @export

simulateTwoCauseFineGrayModel <- function(nobs, beta1, beta2, X = NULL, u.min = 0, u.max, p = 0.5,
                                          returnX = FALSE) {

  #warnings:
  if(length(beta1) != length(beta2)) stop("Dimension of beta1 and beta2 should be the same")
  ncovs <- length(beta1)

  if(is.null(X)) {
    X <- matrix(rnorm(nobs * ncovs), nrow = nobs)
    returnX <- TRUE
  }

  #- Generate indicator for cause
  c.ind <- 1 + rbinom(nobs, 1, prob = (1 - p)^exp(X %*% beta1))

  #Conditional on cause indicators, we simulate the model.
  ftime <- numeric(nobs)

  eta1 <- X[c.ind == 1, ] %*% beta1
  eta2 <- X[c.ind == 2, ] %*% beta2


  u1 <- runif(length(eta1))
  t1 <- -log(1 - (1 - (1 - u1 * (1 - (1 - p)^exp(eta1)))^(1 / exp(eta1))) / p)
  t2 <- rexp(length(eta2), rate = exp(eta2))
  ci <- runif(nobs, min = u.min, max = u.max) #simulate censoring times

  ftime[c.ind == 1] <- t1
  ftime[c.ind == 2] <- t2
  ftime <- pmin(ftime, ci)
  fstatus <- ifelse(ftime == ci, 0, 1)
  fstatus <- fstatus * c.ind

  out         <- list()
  out$ftime   <- ftime
  out$fstatus <- fstatus
  if(returnX)  out$X <- X

  return(out)
}
