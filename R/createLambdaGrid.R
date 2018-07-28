#' Creates Lambda Solution Path for crrBAR
#'
#' @description Creates a grid of \code{lambda} for \code{crrBAR}. This is an internal function.
#'
#' @param ftime A vector of event/censoring times.
#' @param fstatus A vector with unique code for each event type and a separate code for censored observations.
#' @param XX Standardized verion of \code{X} (already standardized in \code{crrBAR})
#' @param uuu Vector of IPCW weights calculated in \code{crrBAR}
#' @param lam.min Numeric: smallest value of lambda if performing grid search
#' @param nlambda Numeric: number of \code{lambda} values if performing grid search  (default is 25)
#' @param log Logical: Whether or not the grid search is log10 spaced (default is \code{TRUE})
#'
#' @details Not intended for use by users. This function is called in \code{crrBAR} if \code{lambda} is missing.
#' @export
createLambdaGrid <- function(ftime, fstatus, XX, uuu,
                        lam.min, nlambda, log = TRUE) {
  n    <- dim(XX)[1]
  p    <- dim(XX)[2]
  eta0 <- rep(0, n)
  sw   <- .C("getScoreAndHessian", as.double(ftime), as.integer(fstatus), as.double(XX),
             as.integer(p), as.integer(n), as.double(uuu), as.double(eta0),
             double(n), double(n), double(1), PACKAGE = "crrBAR")
  score0 <- sw[[8]]
  w0 <- sw[[9]]
  r0 <- ifelse(w0 == 0, 0, score0 / w0)
  z <- eta0 + r0
  l.max <- max(t(w0 * z) %*% XX) / n
  l.min <- lam.min
  if(log) {
    lambda <- 10^(seq(log10(l.max), log10(l.min * l.max), length = nlambda))
  } else {
    lambda <- seq(l.max, l.min * l.max, length = nlambda)
  }
  return(lambda)
}
