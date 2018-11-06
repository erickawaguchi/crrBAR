#include <math.h>
#include <Rmath.h>
#include <string.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R.h>
#include <R_ext/Applic.h>
#include <stdlib.h>
#define LEN sizeof(double)

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
  double grad, hess, shift, likli, s0, si, l1, bhat;
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

      // calculate xwr and xwx & update beta_j
      for (int j = 0; j < p; j++) {
        if(a[j] == 0) b[l * p + j] = 0;

        else {
          //Calculate the jth component of gradient and hessian (with previous values updated)
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
          }

          for ( int j3 = 0; j3 < n; j3++){
            if (w[j3] == 0) r[j3] = 0;
            else r[j3] = st[j3] / w[j3];
          }

          int nj = n * j;
          double grad = 0;
          for (int i = 0; i < n; i++) grad += x[nj + i] * st[i];

          double hess = 0;
          for (int i = 0; i < n; i++) hess += w[i] * pow(x[nj + i], 2);


          l1 = lam[l];
          bhat = a[j];
          //New beta_j update
          b[l * p + j] = newBarL0(hess / n, grad / n, bhat, l1 / n);
        }
        // Update r
        shift = b[l * p + j] - a[j];
        if (shift != 0) {
          for (int i = 0; i < n; i++) {
            si = shift * x[j * n + i]; //low-rank update of beta
            r[i] -= si;
            eta[i] += si; //update eta using new value of beta_j
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
