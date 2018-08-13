#include <math.h>
#include <Rmath.h>
#include <string.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R.h>
#include <R_ext/Applic.h>
#include <stdlib.h>
#define LEN sizeof(double)

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
    st[i]+=1;
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
        zb += b[j1] *x[n * j1 + j];

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
// See Simon et al. 2011 (pg. 4)
double ridge(double z, double l1, double v) {
  // s is the sign function
  double s = 0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  return(s * fabs(z) / (v + l1));
  //return(s * fabs(l2 * z) / (l2 * v + l1));
}

//double lasso(double z, double l1, double l2, double v) {
//    double s=0;
//    if (z > 0) s = 1;
//    else if (z < 0) s = -1;
//    if (fabs(z) <= l1) return(0);
//    else return(s*(fabs(z)-l1)/(v+l2));
//}

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
SEXP cleanupCRR(double *a, double *eta, double *wye, double *st, double *w, SEXP beta, SEXP Dev, SEXP iter, SEXP residuals, SEXP score, SEXP hessian) {
  Free(a);
  Free(eta);
  Free(wye);
  Free(st);
  Free(w);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 6));
  SET_VECTOR_ELT(res, 0, beta); //coefficient estimates
  SET_VECTOR_ELT(res, 1, Dev); //deviance = -2*loglik
  SET_VECTOR_ELT(res, 2, iter); //iterations until convergence
  SET_VECTOR_ELT(res, 3, residuals); //residuals
  SET_VECTOR_ELT(res, 4, score); //gradient
  SET_VECTOR_ELT(res, 5, hessian); //hessian
  UNPROTECT(7);
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
  SEXP res, beta, Dev, iter, residuals, score, hessian;
  PROTECT(beta = allocVector(REALSXP, p));
  double *b = REAL(beta);
  for (int j = 0; j < p; j++) b[j] = 0;
  PROTECT (score = allocVector(REALSXP, n));
  double *s = REAL(score);
  for (int i = 0; i < n; i++) s[i] = 0;
  PROTECT (hessian = allocVector(REALSXP, n));
  double *h = REAL(hessian);
  for (int i = 0; i <  n; i++) h[i] = 0;
  PROTECT(residuals = allocVector(REALSXP, n));
  double *r = REAL(residuals);
  PROTECT(Dev = allocVector(REALSXP, 1));
  for (int i = 0; i < 1; i++) REAL(Dev)[i] = 0;
  PROTECT(iter = allocVector(INTSXP, 1));
  for (int i = 0; i < 1; i++) INTEGER(iter)[i] = 0;
  double *a = Calloc(p, double); // Beta from previous iteration
  for (int j = 0; j < p; j++) a[j] = 0; // Beta0 from previous iteration
  double *st = Calloc(n, double);
  for (int i = 0; i < n; i++) st[i]=0;
  double *w = Calloc(n, double);
  for ( int i = 0; i < n; i++) w[i]=0;
  double *x = REAL(x_);
  double *t2 = REAL(t2_);
  double *wt = REAL(wt_);
  int *ici = INTEGER(ici_);
  double lam = REAL(lambda_)[0];
  double esp = REAL(esp_)[0];
  int max_iter = INTEGER(max_iter_)[0];
  double *m = REAL(multiplier);
  double *eta = Calloc(n, double);
  for (int i = 0; i < n; i++) eta[i] = 0;
  double *wye = Calloc(n, double);
  double xwr, xwx, u, v, l1, shift, likli, s0, si;
  int converged;

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
      xwr = wcrossprod(x, r, w, n, j); // h_j * g_j
      xwx = wsqsum(x, w, n, j); // h_j
      u   = xwr / n + (xwx / n) * a[j]; // z in paper
      v   = xwx / n;

      // Update b_j
      l1 = lam * m[j] / n; //divide by n since we are minimizing the following: -(1/n)l(beta) + lambda * p(beta)

      //Change penalty to ridge regression
      b[j] = ridge(u, l1, v);

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
      s[i] = st[i];
      h[i] = w[i];
    }
    if (converged)  break;
    //for converge
  } //for while loop

  res = cleanupCRR(a, eta, wye, st, w, beta, Dev, iter, residuals, score, hessian);
  return(res);
}

// GROUP BAR
//Penalized grouped ridge
void groupRidge(double *b, double *x, double *r, double *w, double *eta, int g, int *K1, int n, int p, double lam, double *a) {

  // Initialize
  int K = K1[g + 1] - K1[g];
  double *z = Calloc(K, double);
  double *v = Calloc(K, double);

  //Find group norm
  for (int j = K1[g]; j < K1[g + 1]; j++){
    z[j - K1[g]] = wcrossprod(x, r, w, n, j) / n + wsqsum(x, w, n, j) / n * a[j];
    v[j - K1[g]] = wsqsum(x, w, n, j) / n;}
  double z_norm = norm(z, K);

  // Update b
  double len = 0;
  for (int j = K1[g]; j < K1[g + 1]; j++) {
    len = ridge(z_norm, lam, v[j - K1[g]]);
    if ((len != 0) | (a[K1[g]] != 0)) {
      b[p + j] = len * z[j - K1[g]] / z_norm;
      double shift = b[p + j] - a[j];
      for (int i=0; i < n; i++) {
        double si = shift * x[j * n + i];
        r[i] -= si;
        eta[i] += si;
      }//for i=0 to n loop
    }//if loop
  }

  Free(z);
  Free(v);
}


//start group cordinate descent
SEXP ccd_gridge(SEXP x_, SEXP t2_, SEXP ici_, SEXP wt_, SEXP K1_, SEXP penalty_, SEXP lambda_, SEXP esp_, SEXP max_iter_, SEXP multiplier) {
  //Declaration
  int n = length(t2_);
  int p = length(x_) / n;
  int J = length(K1_) - 1; //Number of groups + 1
  int *K1 = INTEGER(K1_);
  double nullDev;
  SEXP res, beta, Dev, iter, residuals, score, hessian;
  PROTECT(beta = allocVector(REALSXP, p));
  double *b = REAL(beta);
  for (int j = 0; j < p; j++) b[j] = 0;
  PROTECT (score = allocVector(REALSXP, n));
  double *s = REAL(score);
  for (int i = 0; i < n; i++) s[i] = 0;
  PROTECT (hessian = allocVector(REALSXP, n));
  double *h = REAL(hessian);
  for (int i = 0; i <  n; i++) h[i] = 0;
  PROTECT(residuals = allocVector(REALSXP, n));
  double *r = REAL(residuals);
  PROTECT(Dev = allocVector(REALSXP, 1));
  for (int i = 0; i < 1; i++) REAL(Dev)[i] = 0;
  PROTECT(iter = allocVector(INTSXP, 1));
  for (int i = 0; i < 1; i++) INTEGER(iter)[i] = 0;
  double *a = Calloc(p, double); // Beta from previous iteration
  for (int j = 0; j < p; j++) a[j] = 0; // Beta0 from previous iteration
  double *st = Calloc(n, double);
  for (int i = 0; i < n; i++) st[i]=0;
  double *w = Calloc(n, double);
  for ( int i = 0; i < n; i++) w[i]=0;
  double *x = REAL(x_);
  double *t2 = REAL(t2_);
  double *wt = REAL(wt_);
  int *ici = INTEGER(ici_);
  double lam = REAL(lambda_)[0];
  double esp = REAL(esp_)[0];
  int max_iter = INTEGER(max_iter_)[0];
  double *m = REAL(multiplier);
  double *eta = Calloc(n, double);
  for (int i = 0; i < n; i++) eta[i] = 0;
  double *wye = Calloc(n, double);
  double l1, likli, s0;
  int converged;

  //end of declaration;

  //initialization
  nullDev = -2 * getLogLikelihood(t2, ici, x, p, n, wt, a); // Calculate null deviance at beta = 0

  for (int j = 0; j < p; j++) a[j] = b[j];

  while (INTEGER(iter)[0] < max_iter) {
    if (REAL(Dev)[0] - nullDev > 0.99 * nullDev) break;

    INTEGER(iter)[0]++;

    // This part is the same as regular BAR
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

    // Here is where we take into account grouping info
    // For each group
    for (int g = 0; g < J; g++) {
      // Update b_j
      l1 = lam * m[g] * sqrt(K1[g + 1] - K1[g]) / n; //divide by n since we are minimizing the following: -(1/n)l(beta) + lambda * p(beta)
      groupRidge(b, x, r, w, eta, g, K1, n, p, l1, a); //use groupRidge penalty defined above
    }

    converged = checkConvergence(b, a, esp, p);
    for (int i = 0; i < p; i++)
      a[i] = b[i];

    //Calculate deviance
    REAL(Dev)[0] = -2 * getLogLikelihood(t2, ici, x, p, n, wt, a);

    for (int i = 0; i < n; i++){
      s[i] = st[i];
      h[i] = w[i];
    }
    if (converged)  break;
    //for converge
  } //for while loop

  res = cleanupCRR(a, eta, wye, st, w, beta, Dev, iter, residuals, score, hessian);
  return(res);
}
