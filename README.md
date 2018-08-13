# crrBAR

Broken Adaptive Ridge (BAR) regression for Competing Risks Regression.

Introduction
============

crrBAR is an `R` package for performing L_0-based regressions for the popular Fine-Gray model for competing risks data.

Dependencies
============
 * `survival`
 

Getting Started
===============
1. On Windows, make sure [RTools](https://CRAN.R-project.org/bin/windows/Rtools/) is installed.
2. In R, use the following commands to download and install crrBAR:

  ```r
  install.packages("devtools")
  library(devtools)
  install_github("erickawaguchi/crrBAR")
  ```

3. To perform L_0-penalized regression, use the following commands in R:
  ```r
  library(crrBAR)
  #Assume cause of interest of fstatus = 1.
  fit <- crrBAR(ftime, fstatus, X, failcode = 1, cencode = 0, lambda = log(sum(fstatus == 1)) / 2, xi = 1 / 2)
  fit$coef #Extract coefficients
  ```
  
Examples
========
 ```r
set.seed(10)
ftime <- rexp(200)
fstatus <- sample(0:2, 200, replace = TRUE)
cov <- matrix(runif(1000), nrow = 200)
dimnames(cov)[[2]] <- c('x1','x2','x3','x4','x5')
fit <- crrBAR(ftime, fstatus, cov, lambda = log(sum(fstatus == 1)) / 2, xi = 1 / 2)
fit$coef
 ```
 
Development
===========
crrBAR is being developed in R Studio. If there are any questions or comments please email me at erickawaguchi[at]ucla.edu.

