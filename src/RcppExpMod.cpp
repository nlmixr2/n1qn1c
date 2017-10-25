#include "evaluate.h"

Rcpp::EvalBase *fev = NULL;                  // pointer to abstract base class
Rcpp::EvalBase *gev = NULL;                  // pointer to abstract base class

typedef void (*S2_fp) (int *, int *, double *, double *, double *, int *, float *, double *);
extern "C" void n1qn1_ (S2_fp simul, int n[], double x[], double f[], double g[], double var[], double eps[],
                        int mode[], int niter[], int nsim[], int imp[], int lp[], double zm[], int izs[], float rzs[], double dzs[]);

unsigned int n1qn1_calls = 0, n1qn1_grads = 0;
int n1qn1_fprint = 0;
static void fwrap(int *ind, int *n, double *x, double *f, double *g, int *ti, float *tr, double *td)
{
  int i;
  Rcpp::NumericVector par(*n), ret(*n);
  for (i = 0; i < *n; i++) par[i] = x[i];
        
  if (*ind==2 || *ind==4) {
    n1qn1_calls++;
    ret = fev->eval(par);
    if (n1qn1_fprint){
      Rprintf("%3d:%#14.8g:", n1qn1_calls, ret[0]);
      for (i = 0; i < *n; i++) Rprintf(" %#8g", x[i]);
      Rprintf("\n");
    }
    *f = ret[0];
  }
  if (*ind==3 || *ind==4) {
    n1qn1_grads++;
    ret = gev->eval(par);
    for (i = 0; i < *n; i++) g[i] = ret[i];
  }
}


RcppExport SEXP
n1qn1_wrap(
           SEXP fSEXP, SEXP gSEXP, SEXP rhoSEXP, SEXP xSEXP, SEXP epsSEXP, 
           SEXP nSEXP, SEXP modeSEXP, SEXP niterSEXP, SEXP nsimSEXP, SEXP impSEXP,
           SEXP nzmSEXP, SEXP zmSEXP, SEXP fprint_sexp) {
  BEGIN_RCPP
    n1qn1_calls=0;
  n1qn1_grads=0;
  n1qn1_fprint = INTEGER(fprint_sexp)[0];
  if (TYPEOF(fSEXP) == EXTPTRSXP){
    fev = new Rcpp::EvalCompiled(fSEXP, rhoSEXP); // xptr
  } else {
    fev = new Rcpp::EvalStandard(fSEXP, rhoSEXP); // Standard evaulation
  }
  if (TYPEOF(gSEXP) == EXTPTRSXP){
    gev = new Rcpp::EvalCompiled(gSEXP, rhoSEXP); // xptr
  } else {
    gev = new Rcpp::EvalStandard(gSEXP, rhoSEXP); // Standard evaulation
  }
  int i, j, k =0;
    
  int n, mode, niter, nsim, imp, lp=6, nzm;
  n = INTEGER(nSEXP)[0];
  mode = INTEGER(modeSEXP)[0];
  niter = INTEGER(niterSEXP)[0];
  nsim = INTEGER(nsimSEXP)[0];
  imp = INTEGER(impSEXP)[0];
  nzm = INTEGER(nzmSEXP)[0];

  double f, eps;
  double *x = new double[n];
  double *g = new double[n];
  double *var = new double[n];
  double *zm = new double[nzm];

  int izs[1]; float rzs[1]; double dzs[1];
  for (i=0; i<n; i++) x[i] = REAL(xSEXP)[i];
  for (i=0; i<nzm; i++) zm[i] = REAL(zmSEXP)[i];
  eps = REAL(epsSEXP)[0];
  for (i=0; i<n; i++) var[i] = .1;

  n1qn1_(fwrap,&n,x,&f,g,var,&eps,
         &mode,&niter,&nsim,&imp,&lp,zm,izs,rzs,dzs);
        
  Rcpp::NumericVector par(n);
  for (i=0; i<n; i++) par[i] = x[i];
  Rcpp::NumericVector hess(nzm);
  // On input this is hessian
  // On output this is H = LDL'
  // Triangular matrix is paramterized by column instead of row.
  using namespace arma;
  mat L = mat(n,n);
  mat D = mat(n,n);
  mat H = mat(n,n);
  // NumericVector zms(nzm);
  // for (i = 0; i < nzm; i++) zms[i]=zm[i];
  L.zeros();
  D.zeros();
  k =0;
  for (i=0; i<n; i++){
    for (j=i; j<n; j++){
      if (i == j){
        D(i,i)=zm[k];
        L(i,i)=1;
      } else {
        L(j,i)=zm[k];
      }
      k++;
    }
  }
  H = L*D*L.t();
  k = 0;
  for (i=0; i<n; i++){
    for (j=i; j<n; j++){
      hess[k]=H(j,i);
      k++;
    }
  }

  delete[] x;
  delete[] g;
  delete[] var;
  delete[] zm;

  return Rcpp::List::create(Rcpp::Named("value") = f,
                            Rcpp::Named("par") = par,
                            // Rcpp::Named("L") = L,
                            // Rcpp::Named("D") = D,
                            Rcpp::Named("H") = H,
                            // Rcpp::Named("zm")=zms,
                            Rcpp::Named("c.hess") = hess,
			    Rcpp::Named("n.fn") = n1qn1_calls,
			    Rcpp::Named("n.gr") = n1qn1_grads);
  END_RCPP
 }
