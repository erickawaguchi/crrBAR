#include <math.h>
#include <Rmath.h>
#include <string.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R.h>
#include <R_ext/Applic.h>
#include <stdlib.h>
#define LEN sizeof(double)

double sgn(double z) {
  double s = 0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  return(s);
}

void getScoreAndHessian(double *t2, int *ici, double *x, int *ncov, int *nin, double *wt, double *eta, double *st, double *w, double *lik)
{
  const int np = ncov[0], n = nin[0];
  int i, j2, k, j;
  double likli = 0, s0, *a, wye[n];
  // a is the matrix of covariates
  //s1 is the score
  a = (double*)malloc(n * np * LEN);
  //initialization

  for (int i = 0; i < n; i++)
    for (int k = 0; k < np; k++)
      *(a + i * np + k) = x[i + n * k];

  //start here
  for (i = 0; i < n; i++)
  {
    if (ici[i] != 1) continue;
    likli += eta[i];
    st[i] += 1;
    // score
    s0 = 0;
    for (j = 0; j < n; j++)
      wye[j] = 0;
    for (k = 0; k < n; k++)
    {

      if (t2[k] < t2[i] && ici[k] <= 1) continue;//leave those out of risk set
      if (t2[k] >= t2[i])
        wye[k] = exp(eta[k]);
      else
        wye[k] = exp(eta[k]) * wt[i] / wt[k];

      s0 += wye[k];
    }
    for (j2 = 0; j2 < n; j2++){
      st[j2] += -wye[j2] / s0;
      w[j2] += wye[j2] / s0 - pow(wye[j2], 2) / (s0 * s0);
    }
    likli -= log(s0);
    //end of score
  }

  *lik = likli;
  free(a);
}

//Standardize design matrix
SEXP standardize(SEXP X_) {
  // Declarations
  int n = nrows(X_);
  int p = ncols(X_);
  SEXP XX_, c_, s_;
  PROTECT(XX_ = allocMatrix(REALSXP, n, p));
  PROTECT(c_ = allocVector(REALSXP, p));
  PROTECT(s_ = allocVector(REALSXP, p));
  double *X = REAL(X_);
  double *XX = REAL(XX_);
  double *c = REAL(c_);
  double *s = REAL(s_);

  for (int j = 0; j < p; j++) {

    // Center (Calculate mean and subtract)
    c[j] = 0;
    for (int i = 0; i < n; i++) {
      c[j] += X[j * n + i];
    }
    c[j] = c[j] / n;
    for (int i = 0; i < n; i++) XX[j * n + i] = X[j * n + i] - c[j];

    // Scale (Calculate sdev and divide)
    s[j] = 0;
    for (int i = 0; i < n; i++) {
      s[j] += pow(XX[j * n + i], 2);
    }
    s[j] = sqrt(s[j] / n);
    for (int i = 0; i < n; i++) XX[j * n + i] = XX[j * n + i] / s[j];
  }

  // Return list
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(res, 0, XX_); // Standardized design matrix
  SET_VECTOR_ELT(res, 1, c_); // Mean
  SET_VECTOR_ELT(res, 2, s_); // Standard deviations
  UNPROTECT(4);
  return(res);
}

//From penalty file

double getLogLikelihood(double *t2, int *ici, double *x, int ncov, int nin, double *wt, double *b)
{
  // Look at Eq (2) from Fu et al. 2017.
  int i,j, j1;
  const int p = ncov,  n = nin;
  double likli = 0, zb, s0;

  for (i = 0; i < n; i++)
  {
    if (ici[i] != 1) continue;
    //first part of log-partial likelihood
    for (j = 0; j < p; j++)
      likli += b[j] * x[n * j + i];

    //second part of log-partial likelihood
    s0 = 0;
    for (j = 0; j < n; j++)
    {
      if (t2[j] < t2[i] && ici[j] <= 1) continue;
      zb = 0.0;

      for (j1 = 0; j1 < p; j1 ++)
        zb += b[j1] * x[n * j1 + j];

      if (t2[j] >= t2[i])
        s0 += exp(zb);
      else
        s0 += exp(zb) * wt[i] / wt[j];
    }

    likli -= log(s0);
  }
  return likli;
}


SEXP evalLogLikelihood(SEXP x_, SEXP t2_, SEXP ici_, SEXP wt_, SEXP beta_) {
  //Declaration
  int n = length(t2_);
  int p = length(x_) / n;
  //int L = length(lambda);

  // initialize
  double *x = REAL(x_);
  double *t2 = REAL(t2_);
  double *wt = REAL(wt_);
  double *b = REAL(beta_);
  int *ici = INTEGER(ici_);
  double loglik = getLogLikelihood(t2, ici, x, p, n, wt, b);
  return(ScalarReal(loglik));
}


// Criterion for convergence: All coefficients must pass the following |(b_new - b_old) / b_old| < eps
int checkConvergence(double *beta, double *beta_old, double eps, int p) {
  int converged = 1;
  for (int j = 0; j < p; j++) {
    if (fabs((beta[j] - beta_old[j]) / beta_old[j]) > eps) {
      converged = 0;
      break;
    }
  }
  return(converged);
}

// Create penalty function for ridge regression. Lasso is here as outline
// Edit: Added into BAR function directly.
// See Simon et al. 2011 (pg. 4)
double ridge(double z, double l1, double v) {
  return(z / (v + l1));
}

// Euclidean norm (||x||_2 = \sqrt{\sum x_i^2})
double norm(double *x, int p) {
  double x_norm = 0;
  for (int j = 0; j < p; j++) x_norm = x_norm + pow(x[j], 2);
  x_norm = sqrt(x_norm);
  return(x_norm);
}

// Weighted cross product of y with jth column of x
double wcrossprod(double *X, double *y, double *w, int n, int j) {
  int nn = n * j;
  double val = 0;
  for (int i = 0; i < n; i++) val += X[nn + i] * y[i] * w[i];
  return(val);
}

// Weighted sum of squares of jth column of X
double wsqsum(double *X, double *w, int n, int j) {
  int nn = n * j;
  double val = 0;
  for (int i = 0; i < n; i++) val += w[i] * pow(X[nn + i], 2);
  return(val);
}


///////////////////////////////////////////////////////////////////////////////////////
SEXP cleanupCRR(double *a, double *eta, double *wye, double *st, double *w, double *diffBeta,
                SEXP beta, SEXP Dev, SEXP iter, SEXP residuals, SEXP score, SEXP hessian, SEXP linpred) {
  Free(a);
  Free(eta);
  Free(wye);
  Free(st);
  Free(w);
  Free(diffBeta);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 7));
  SET_VECTOR_ELT(res, 0, beta); //coefficient estimates
  SET_VECTOR_ELT(res, 1, Dev); //deviance = -2*loglik
  SET_VECTOR_ELT(res, 2, iter); //iterations until convergence
  SET_VECTOR_ELT(res, 3, residuals); //residuals
  SET_VECTOR_ELT(res, 4, score); //gradient
  SET_VECTOR_ELT(res, 5, hessian); //hessian
  SET_VECTOR_ELT(res, 6, linpred); //hessian
  UNPROTECT(8);
  return(res);
}


//////////////////////////////////////////////////////////////////////////////////
//start cordinate descent

//x_ = design matrix
//t2_ = failtime
//ici_ = censoring vector
//wt_ = weight (uuu)
//lambda = tuning parameter
//esp_ = epsilon (thershold)
//max_iter_ = max iterations
//multiplier = penalty.factor (what to reweight lambda by)

SEXP ccd_ridge(SEXP x_, SEXP t2_, SEXP ici_, SEXP wt_, SEXP lambda_,
               SEXP esp_, SEXP max_iter_, SEXP multiplier) {

  //Declaration
  int n = length(t2_);
  int p = length(x_) / n;
  double nullDev;

  //Output
  SEXP res, beta, Dev, iter, residuals, score, hessian, linpred;
  PROTECT(beta = allocVector(REALSXP, p));
  double *b = REAL(beta);
  for (int j = 0; j < p; j++) b[j] = 0;
  PROTECT (score = allocVector(REALSXP, n));
  double *s = REAL(score);
  for (int i = 0; i < n; i++) s[i] = 0;
  PROTECT (hessian = allocVector(REALSXP, n));
  double *h = REAL(hessian);
  for (int i = 0; i < n; i++) h[i] = 0;

  PROTECT(residuals = allocVector(REALSXP, n));
  double *r = REAL(residuals);
  PROTECT(Dev = allocVector(REALSXP, 1));
  for (int i = 0; i < 1; i++) REAL(Dev)[i] = 0;
  PROTECT(iter = allocVector(INTSXP, 1));
  for (int i = 0; i < 1; i++) INTEGER(iter)[i] = 0;
  PROTECT(linpred = allocVector(REALSXP, n));
  double *lp = REAL(linpred);
  for (int i = 0; i <  n; i++) lp[i] = 0;

  //Intermediate quantities for internal use (must be freed afterwards!)
  double *a = Calloc(p, double); // Beta from previous iteration
  for (int j = 0; j < p; j++) a[j] = 0;
  double *st = Calloc(n, double);
  for (int i = 0; i < n; i++) st[i] = 0;
  double *w = Calloc(n, double);
  for ( int i = 0; i < n; i++) w[i] = 0;
  double *eta = Calloc(n, double);
  for (int i = 0; i < n; i++) eta[i] = 0;
  double *diffBeta = Calloc(p, double);
  for (int j = 0; j < p; j++) diffBeta[j] = 1;
  double *wye = Calloc(n, double);
  double xwr, xwx, l1, shift, likli, s0, si, delta;
  int converged;

  //Pointers
  double *x = REAL(x_);
  double *t2 = REAL(t2_);
  double *wt = REAL(wt_);
  int *ici = INTEGER(ici_);
  double lam = REAL(lambda_)[0];
  double esp = REAL(esp_)[0];
  int max_iter = INTEGER(max_iter_)[0];
  double *m = REAL(multiplier);


  //end of declaration;

  //initialization
  nullDev = -2 * getLogLikelihood(t2, ici, x, p, n, wt, a); // Calculate null deviance at beta = 0

  for (int j = 0; j < p; j++) a[j] = b[j];

  while (INTEGER(iter)[0] < max_iter) {
    if (REAL(Dev)[0] - nullDev > 0.99 * nullDev) break;

    INTEGER(iter)[0]++;

    //calculate score and w
    likli = 0;
    for (int i = 0; i < n; i++){
      st[i] = 0;
      w[i]= 0;
    }
    for (int i = 0; i < n; i++)
    {
      if (ici[i] != 1) continue;
      likli += eta[i];
      st[i] += 1;

      // score
      s0 = 0;
      for (int j1 = 0; j1 < n; j1++)
        wye[j1] = 0;
      for (int k = 0; k < n; k++)
      {
        if (t2[k] < t2[i] && ici[k] <= 1) continue;

        if (t2[k] >= t2[i])
          wye[k] = exp(eta[k]);
        else
          wye[k] = exp(eta[k]) * wt[i] / wt[k];
        s0 += wye[k];
      }
      for (int j2 = 0; j2 < n; j2++){
        st[j2] += -wye[j2] / s0;
        w[j2] += wye[j2] / s0 - pow(wye[j2], 2) / (s0 * s0);
      }
      likli -= log(s0);
    }

    for ( int j3 = 0; j3 < n; j3++){
      if (w[j3] == 0) r[j3] = 0;
      else r[j3] = st[j3] / w[j3];
    }

    // calculate xwr and xwx & update beta_j
    for (int j = 0; j < p; j++) {
      xwr = -wcrossprod(x, r, w, n, j); // jth component of gradient [l'(b)]
      xwx = wsqsum(x, w, n, j); // jth component of hessian [l''(b)]
      l1 = lam * m[j] / n; //divide by n since we are minimizing the following: -(1/n)l(beta) + lambda * p(beta)
      delta =  -(xwr / n + a[j] * l1) / (xwx / n + l1);

      //u   = xwr / n + (xwx / n) * a[j]; // z in paper
      //v   = xwx / n;

      // Update b_j

      // Do one dimensional ridge update.
      // Employ trust region as in Genkin et al. (2007) for quadratic approximation.
      b[j] = a[j] + sgn(delta) * fmin(fabs(delta), diffBeta[j]);
      diffBeta[j] = fmax(2 * fabs(delta), diffBeta[j] / 2);

      // Update r
      shift = b[j] - a[j];
      if (shift != 0) {
        for (int i = 0; i < n; i++) {
          si = shift * x[j * n + i];
          r[i] -= si;
          eta[i] += si;
        }
      } //end shift

    } //for j = 0 to (p - 1)
    // Check for convergence
    converged = checkConvergence(b, a, esp, p);
    for (int i = 0; i < p; i++)
      a[i] = b[i];

    //Calculate deviance
    REAL(Dev)[0] = -2 * getLogLikelihood(t2, ici, x, p, n, wt, a);

    for (int i = 0; i < n; i++){
      lp[i] = eta[i];
      s[i] = st[i];
      h[i] = w[i];
    }
    if (converged)  break;
    //for converge
  } //for while loop

  res = cleanupCRR(a, eta, wye, st, w, diffBeta, beta, Dev, iter, residuals, score, hessian, linpred);
  return(res);
}

////////////////////////////////////////////////////
// FAST L_0-BAR
////////////////////////////////////////////////////
int checkFastBarConvergence(double *beta, double *beta_old, double eps, int l, int p) {
  int converged = 1;
  for (int j = 0; j < p; j++) {
    if (fabs((beta[l * p + j] - beta_old[j])) > eps) {
      converged = 0;
      break;
    }
  }
  return(converged);
}

SEXP cleanupNewCRR(double *a, double *e,  double *eta, double *wye, double *st, double *w,
                   SEXP beta, SEXP Dev, SEXP iter, SEXP residuals, SEXP score, SEXP hessian, SEXP linpred, SEXP converged) {
  // Free up all intermediate step variables
  Free(a);
  Free(e);
  Free(eta);
  Free(wye);
  Free(st);
  Free(w);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 8));
  SET_VECTOR_ELT(res, 0, beta); //coefficient estimates
  SET_VECTOR_ELT(res, 1, Dev); //deviance = -2*loglik
  SET_VECTOR_ELT(res, 2, iter); //iterations until convergence
  SET_VECTOR_ELT(res, 3, residuals); //residuals
  SET_VECTOR_ELT(res, 4, score); //gradient
  SET_VECTOR_ELT(res, 5, hessian); //hessian
  SET_VECTOR_ELT(res, 6, linpred); //hessian
  SET_VECTOR_ELT(res, 7, converged); //check convergence
  UNPROTECT(9);
  return(res);
}

double newBarL0(double h, double g, double b, double l) {
  double tmp;
  double s = 0;
  tmp = h * b + g;
  if (tmp > 0) s = 1;
  else if (tmp < 0) s = -1;
  if (fabs(tmp) < 2 * sqrt(h * l)) return(0);
  else return((tmp + s * sqrt(pow(tmp, 2) - 4 * l * h)) / (2 * h));
}

//double newBarL1(double h, double g, double b, double l) {
//  double tmp;
//  tmp = h * b + g;
//  if (fabs(tmp) < 2 * sqrt(l * h)) return(0);
//  else return((tmp + sqrt(pow(tmp, 2) - 4 * l * h)) / 2 * h);
//}

// CCD BAR L_0 or L_1
SEXP ccd_bar(SEXP x_, SEXP t2_, SEXP ici_, SEXP wt_, SEXP lambda_,
                 SEXP esp_, SEXP max_iter_, SEXP beta0_) {

  //Declaration
  int n = length(t2_);
  int p = length(x_) / n;
  int L = length(lambda_);
  double nullDev;

  //Output
  SEXP res, beta, Dev, iter, residuals, score, hessian, converged, linpred;
  PROTECT(beta = allocVector(REALSXP, L * p));
  double *b = REAL(beta);
  for (int j = 0; j < (L * p); j++) b[j] = 0;
  PROTECT (score = allocVector(REALSXP, L * n));
  double *s = REAL(score);
  for (int i = 0; i < (L * n); i++) s[i] = 0;
  PROTECT (hessian = allocVector(REALSXP, L * n));
  double *h = REAL(hessian);
  for (int i = 0; i <  (L * n); i++) h[i] = 0;
  PROTECT(residuals = allocVector(REALSXP, n));
  double *r = REAL(residuals);
  PROTECT(Dev = allocVector(REALSXP, L + 1));
  for (int i = 0; i < (L + 1); i++) REAL(Dev)[i] = 0;
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i = 0; i < L; i++) INTEGER(iter)[i] = 0;
  PROTECT(converged = allocVector(INTSXP, L));
  for (int i = 0; i < L; i++) INTEGER(converged)[i] = 0;
  PROTECT(linpred = allocVector(REALSXP, n));
  double *lp = REAL(linpred);
  for (int i = 0; i <  n; i++) lp[i] = 0;

  //Intermediate quantities for internal use (free this memory)
  double *a = Calloc(p, double); // Beta from previous iteration
  for (int j = 0; j < p; j++) a[j] = 0;
  double *st = Calloc(n, double);
  for (int i = 0; i < n; i++) st[i]=0;
  double *w = Calloc(n, double);
  for ( int i = 0; i < n; i++) w[i]=0;
  double *eta = Calloc(n, double);
  for (int i = 0; i < n; i++) eta[i] = 0;
  double *e = Calloc(n, double);
  for (int i = 0; i < n; i++) e[i] = 0;
  double *wye = Calloc(n, double);
  double xwr, xwx, shift, likli, s0, si, l1, bhat;
  //int converged;

  //Pointers
  double *x = REAL(x_);
  double *t2 = REAL(t2_);
  double *wt = REAL(wt_);
  int *ici = INTEGER(ici_);
  double *lam = REAL(lambda_);
  double esp = REAL(esp_)[0];
  int max_iter = INTEGER(max_iter_)[0];
  double *m = REAL(beta0_);
  //double *e = REAL(eta_);
  //const char *penalty = CHAR(STRING_ELT(penalty_, 0));
  //end of declaration;

  //initialization

  //Initialize beta_0 = beta^(0)

  //Initialize eta with beta_0
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      e[i] += m[j] * x[j * n + i];
    }
  }

  nullDev = -2 * getLogLikelihood(t2, ici, x, p, n, wt, a); // Calculate null deviance at beta = 0
  REAL(Dev)[0] = nullDev; //Store initial loglikelihood as first element in Dev

  // Initialize a and eta using m and e. (Warm starts will be used afterwards)
  for (int j = 0; j < p; j++) a[j] = m[j];
  for (int i = 0; i < n; i++) eta[i] = e[i];

  //Outer loop for each lambda
  for(int l = 0; l < L; l++) {

    //Start algorithm here
    while (INTEGER(iter)[l] < max_iter) {

      if (REAL(Dev)[l + 1] - nullDev > 0.99 * nullDev) break;

      INTEGER(iter)[l]++;

      //calculate score and w
      likli = 0;
      for (int i = 0; i < n; i++){
        st[i] = 0;
        w[i]= 0;
      }
      for (int i = 0; i < n; i++)
      {
        if (ici[i] != 1) continue;
        likli += eta[i];
        st[i] += 1;

        // score
        s0 = 0;
        for (int j1 = 0; j1 < n; j1++)
          wye[j1] = 0;
        for (int k = 0; k < n; k++)
        {
          if (t2[k] < t2[i] && ici[k] <= 1) continue;

          if (t2[k] >= t2[i])
            wye[k] = exp(eta[k]);
          else
            wye[k] = exp(eta[k]) * wt[i] / wt[k];
          s0 += wye[k];
        }
        for (int j2 = 0; j2 < n; j2++){
          st[j2] += -wye[j2] / s0;
          w[j2] += wye[j2] / s0 - pow(wye[j2], 2) / (s0 * s0);
        }
        likli -= log(s0);
      }

      for ( int j3 = 0; j3 < n; j3++){
        if (w[j3] == 0) r[j3] = 0;
        else r[j3] = st[j3] / w[j3];
      }

      // calculate xwr and xwx & update beta_j
      for (int j = 0; j < p; j++) {

        xwr = wcrossprod(x, r, w, n, j); //  g_j
        xwx = wsqsum(x, w, n, j); // h_j
        l1 = lam[l];
        bhat = a[j];
        //New beta_j update
        b[l * p + j] = newBarL0(xwx / n, xwr / n, bhat, l1 / n);

        // Update r
        shift = b[l * p + j] - a[j];
        if (shift != 0) {
          for (int i = 0; i < n; i++) {
            si = shift * x[j * n + i]; //low-rank update of beta
            r[i] -= si;
            eta[i] += si;
          }
        } //end shift

      } // End cyclic coordinate-wise optimization

      // Check for convergence (b = current est. a = old estimate)
      INTEGER(converged)[l] = checkFastBarConvergence(b, a, esp, l, p);

      for (int j = 0; j < p; j++)
        a[j] = b[l * p + j]; //redefine old est as new est

      //Calculate deviance
      REAL(Dev)[l + 1] = -2 * getLogLikelihood(t2, ici, x, p, n, wt, a);
      for (int i = 0; i < n; i++){
        s[l * n + i] = st[i];
        h[l * n + i] = w[i];
        lp[i] = eta[i];
      }
      if (INTEGER(converged)[l])  break;
    } //BAR iterate for l^th lambda
  } // Cycle through all lambda

  res = cleanupNewCRR(a, e, eta, wye, st, w, beta, Dev, iter, residuals, score, hessian, linpred, converged);
  return(res);
}
