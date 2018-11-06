crrBARL0old <- function(ftime, fstatus, X, failcode = 1, cencode = 0,
                      lambda = 0, xi = 0, delta = 0,
                      eps = 1E-6,
                      lam.min = ifelse(dim(X)[1] > dim(X)[2], 0.001, 0.05),
                      nlambda = 25,
                      log = TRUE,
                      max.iter = 1000){

  ## Error checking
  if(xi < 0) stop("xi must be a non-negative number.")
  if(max.iter < 1) stop("max.iter must be positive integer.")
  if(eps <= 0) stop("eps must be a positive number.")
  # Sort time
  n <- length(ftime)
  p <- ncol(X)
  d <- data.frame(ftime = ftime, fstatus = fstatus)
  if (!missing(X)) d$X <- as.matrix(X)
  d        <- d[order(d$ftime), ]
  ftime    <- d$ftime
  cenind   <- ifelse(d$fstatus == cencode, 1, 0)
  fstatus  <- ifelse(d$fstatus == failcode, 1, 2 * (1 - cenind))
  X <- d$X
  u <- do.call('survfit', list(formula = Surv(ftime, cenind) ~ 1,
                               data = data.frame(ftime, cenind)))

  # uuu is weight function (IPCW)
  u <- approx(c(0, u$time, max(u$time) * (1 + 10 * .Machine$double.eps)), c(1, u$surv, 0),
              xout = ftime * (1 - 100 * .Machine$double.eps), method = 'constant',
              f = 0, rule = 2)
  uuu <- u$y

  # Standardize design matrix here
  std    <- .Call("standardize", X, PACKAGE = "crrBAR")
  XX     <- std[[1]]
  center <- std[[2]]
  scale  <- std[[3]]
  nz <- which(scale > 1e-6)
  if (length(nz) != ncol(XX)) XX <- XX[ , nz, drop = FALSE]

  # If lambda is MISSING, create lambda path (be sure lambda is sorted)
  if(missing(lambda)) {
    lambda <- createLambdaGrid(ftime, fstatus, XX, uuu, lam.min = lam.min,
                               nlambda = nlambda, log = log)
  } else {
    if(is.unsorted(lambda)) sort(lambda)
  }
  if(min(lambda) < 0) stop("lambda must be a non-negative number.")
  nlam <- length(lambda)

  ## Reorganize ftime, fstatus, and XX in decreasing order:
  #XX      <- XX[order(ftime, decreasing = TRUE), ]
  #ftime   <- rev(ftime)
  #fstatus <- rev(fstatus)
  #uuu     <- rev(uuu)

  ## Fit the PSH Ridge Model here w/ tuning parameter xi
  ridgeFit   <- .Call("ccd_ridge", XX, as.numeric(ftime), as.integer(fstatus), uuu,
                      xi, eps, as.integer(max.iter),
                      penalty.factor = rep(1, p), d = 2, PACKAGE = "crrBAR")
  ridgeCoef  <- ridgeFit[[1]] / scale #Divide coeff estimates by sdev
  ridgeIter  <- ridgeFit[[3]]

  #Enter BAR Fit here (for different lambdas) keep everything in terms of standardized coefficients
  btmp <- ridgeFit[[1]]
  barFit <- .Call("ccd_bar_alt", XX, as.numeric(ftime), as.integer(fstatus), uuu,
                  as.vector(lambda), as.double(eps), as.integer(max.iter),
                  as.vector(btmp), PACKAGE = "crrBAR")

  #Store results here
  coefMatrix  <- as.matrix(matrix(barFit[[1]], p)) / scale
  logLik      <- as.double(barFit[[2]][-1] / -2)
  logLik.null <- as.double(barFit[[2]][1] / -2)
  iter        <- as.integer(barFit[[3]])
  conv        <- as.integer(barFit[[8]])
  grad        <- as.matrix(barFit[[5]], n)
  hess        <- as.matrix(barFit[[6]], n)

  ## Output
  colnames(coefMatrix) <- round(lambda, 3)
  val <- structure(list(coef = coefMatrix,
                        logLik = logLik,
                        logLik.null = logLik.null,
                        grad = grad,
                        hess = hess,
                        iter = iter,
                        lambda = lambda,
                        converged = conv,
                        ridgeCoef = ridgeCoef,
                        #ridgeIter = ridgeIter,
                        xi = xi,
                        call = sys.call()),
                   class = "crrBAR")

  val
}
