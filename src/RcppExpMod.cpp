#include "evaluate.h"
#include <limits>
using namespace arma;

Rcpp::EvalBase *fev = NULL;                  // pointer to abstract base class
Rcpp::EvalBase *gev = NULL;                  // pointer to abstract base class

typedef void (*S2_fp) (int *, int *, double *, double *, double *, int *, float *, double *);
extern "C" void n1qn1_ (S2_fp simul, int n[], double x[], double f[], double g[], double var[], double eps[],
                        int mode[], int niter[], int nsim[], int imp[], double zm[], int izs[], float rzs[], double dzs[]);

unsigned int n1qn1_calls = 0, n1qn1_grads = 0;
int n1qn1_fprint = 0;
static void fwrap(int *ind, int *n, double *x, double *f, double *g, int *ti, float *tr, double *td)
{
  int i;
  Rcpp::NumericVector par(*n), ret(*n);
  std::copy(&x[0], &x[0]+*n, &par[0]);
  
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
    std::copy(&ret[0],&ret[0]+*n,&g[0]);
    // for (i = 0; i < *n; i++) g[i] = ret[i];
  }
}

uvec lowerTri(mat H, bool diag = false){
  unsigned int d = H.n_rows;
  mat o(d, d, fill::ones);
  if (!diag){
    return find(trimatl(o,-1));
  } else {
    return find(trimatl(o));
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
  std::copy(&(REAL(xSEXP)[0]),&(REAL(xSEXP)[0])+n, &x[0]);
  std::copy(&(REAL(zmSEXP)[0]),&(REAL(zmSEXP)[0])+nzm, &zm[0]);
  eps = REAL(epsSEXP)[0];
  std::fill(&var[0], &var[0]+n, 0.1);
  
  n1qn1_(fwrap,&n,x,&f,g,var,&eps,
         &mode,&niter,&nsim,&imp,zm,izs,rzs,dzs);
        
  Rcpp::NumericVector par(n);
  std::copy(&x[0],&x[0]+n,&par[0]);
  // for (i=0; i<n; i++) par[i] = x[i];
  Rcpp::NumericVector hess(nzm);
  // On input this is hessian
  // On output this is H = LDL'
  // Triangular matrix is paramterized by column instead of row.
  mat L = eye(n,n);
  mat D = mat(n,n,fill::zeros);
  mat H = mat(n,n);
  vec zmV(n*(n+1)/2);
  std::copy(&zm[0], &zm[0]+n*(n+1)/2, zmV.begin());
  H.elem(lowerTri(H,true)) = zmV;
  if (n == 1) H = D;
  else{
    L.elem(lowerTri(H,false)) = H.elem(lowerTri(H,0));
    D.diag() = H.diag();
    H = L*D*L.t();
  }
  // Hessian -> c.hess
  vec hessV = H.elem(lowerTri(H,true));
  std::copy(hessV.begin(),hessV.end(),hess.begin());
  
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

extern "C" void qnbd_(int* indqn, S2_fp, int* n, double* x, double* f, double* g, int* iprint, double* zero, int* napmax, 
		      int* itmax, double* epsf, double* epsg, double* epsx, double* df0, 
		      double* binf, double* binsup, int* nfac, double* trav, int* ntrav, int* itrav, int* nitrav, 
		      int* izs, float* rzs, double* dzs);
/*
  qnbd header translated to english by google translate:

  code for minimizing a regular function under constraints of terminals, 
  to standards modulopt method principle of the algorithm: quasi-newton + projection

  details in the report inria n. 242.1983
  this version allows to test several variants of the algorithm
  by modifying some internal parameters (cf how in the code)
  memory size required of the order of n ** 2/2
  for large problems the code gcbd is better adapted under.

  Programs called:
  - zqnbd effective optimizer projection
  - calmaj update of the hessian
  - update update of choleski factors
  - rlbd, saturated linear search with call list terminals

  subroutine qnbd(indqn,simul,n,x,f,g,iprint,zero,
  &     napmax,itmax,epsf,epsg,epsx,df0,binf,bsup,nfac,
  & trav,ntrav,itrav,nitrav,izs,rzs,dzs)


  indqn indicator of qnbd 
  in entry = 1 standard
  = 2 dh and indic initialises at the beginning of work and itrav
  ifac, f, g initialises  output

  on exit
  if <0 unable to calculate a point better than the starting point
  if = 0 stop requests by the user
  if> 0 we provide a better point than the starting point

  <-10 unsuitable input parameters
  = -6 stop when calculating the descent direction and iter = 1
  = -5 stop when calculating the approximation of the hessian iter = 1
  = -3 Simul anomaly: negative sign at a point or
       f and g were previously calculated
  = -2 failure of the linear search at the first iteration
  = -1 f not defined at the initial point
  = 1 stop on epsg
  = 2 epsf
  = 3 epsx
  = 4 napmax
  = 5 itmax
  = 6 slope in the opposite direction to the gradient too small
  = 7 stop when calculating the descent direction
  = 8 stop when calculating the Hessian approximation
  = 10 stop by failure of the linear search, cause not specified
  = 11 ditto with indsim <0
  = 12 a step too small close to a step too big
  this can result from an error in the gradient
  = 13 too many calls in a linear search
     
  simul see modulopt standards
  n dim of x e
  binf, bsup terminals inf, sup, de dim n 
  x variables to optimize (control) 
  f value of the criterion(s)
  g gradient of fs
  zero close zero machine
  napmax maximum number of simulated calls
  itmax maximum number of descents
  itrav vect work dim nitrav = 2n, decomposes into indic and izig
  nfac number of factorized variables (e if indqn = 2)
  iprint print factor 
  varies from 0 (no impressions) to 3 (many impressions)
  epsx vect dim n precision on x e
  epsf criteria stop on e
  epsg stop if sup a norm2 (g +) / n e
  work vect work dim ntrav
  you need ntrav> n (n + 1) / 2 + 6n
  df0> 0 decay f expected (take 1. default)
  izs, rzs, dzs: cf modulopt standards

  indications on the internal variables a qnbd and zqnbd
  izig serves for storing constraints (active if izag> 1)
  if i do not change each one remove one izig (positive)
  otherwise we add izag
  factorization only if izig is null
  dh hessian estimation dim n (n + 1) / 2 ordered in three pieces
  indic (i) new index of the index i
  indic vect dim n storage order of indices es
  not necessary to initialize if indqn = 1

  Parameters of linear search
  amd, amf param. of the wolfe test. (.7 .1)
  napm max number of calls in the rl (= 15)
*/

RcppExport SEXP
qnbd_wrap(SEXP fSEXP, SEXP gSEXP, SEXP rhoSEXP, SEXP xSEXP, 
	  SEXP lowerSEXP, SEXP upperSEXP, SEXP zeroSEXP, SEXP napmaxSEXP,
	  SEXP itmaxSEXP, SEXP epsfSEXP, SEXP epsgSEXP, SEXP epsxSEXP, SEXP nSEXP,
          SEXP fprint_sexp) {
  BEGIN_RCPP
    n1qn1_calls=0;
  n1qn1_grads=0;
  n1qn1_fprint = INTEGER(fprint_sexp)[0];
    int indqn[1];
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
  indqn[0]=1; // 2 is a restart, not really supported.
  int n;
  n = INTEGER(nSEXP)[0];
  int nitrav = 2*n;
  int ntrav = n * (n + 1) / 2 + 6*n+1;
  double *x = new double[n];
  double *g = new double[n];
  double *lowerIn = REAL(lowerSEXP);
  double *upperIn = REAL(upperSEXP);
  double *lower = new double[n];
  double *upper = new double[n];
  for (int j = n; j--;){
    if (R_FINITE(lowerIn[j])){
      lower[j] = lowerIn[j];
    } else {
      lower[j] = std::numeric_limits<double>::lowest();
    }
    if (R_FINITE(upperIn[j])){
      upper[j] = upperIn[j];
    } else {
      upper[j] = std::numeric_limits<double>::max();
    }
  }
  int *itrav = new int[nitrav];
  double *trav = new double[ntrav];
  double *epsx = new double[n];
  int iprint = 0;
  double zero = REAL(zeroSEXP)[0];
  // maximum number of function calls
  int napmax = INTEGER(napmaxSEXP)[0];
  int itmax = INTEGER(itmaxSEXP)[0];
  double epsf = REAL(epsfSEXP)[0];
  double epsg = REAL(epsgSEXP)[0];
  double epsx0 = REAL(epsxSEXP)[0];
  double *xS = REAL(xSEXP);
  std::copy(&xS[0],&xS[0]+n,&x[0]);
  double df0 = 1.0;
  int nfac = 0;
  double f;
  int izs[1]; float rzs[1]; double dzs[1];
  std::fill(&epsx[0], &epsx[0]+n, epsx0);
  qnbd_(indqn, fwrap, &n,x, &f, g, &iprint, &zero,
	&napmax, &itmax, &epsf, &epsg, epsx, &df0,
	lower, upper, &nfac, trav, &ntrav, itrav,
	&nitrav, izs, rzs, dzs);
  Rcpp::NumericVector par(n);
  std::copy(&x[0],&x[0]+n,&par[0]);
  Rcpp::NumericVector travS(ntrav);
  std::copy(&trav[0],&trav[0]+ntrav,&travS[0]);
  
  Rcpp::CharacterVector ret(1);
  switch (indqn[0]){
  case -6: ret[0] = "Stopped when calculating the descent direction on first iteration."; break;
  case -5: ret[0] = "Stopped when calculating the Hessain approximation on first iteration."; break;
  case -3: ret[0] = "function call anomaly: negative sign at a point or f and g were previously calculated"; break;
  case -2: ret[0] = "Failure of the linear search at the first iteration"; break;
  case -1: ret[0] = "f not defined at initial point"; break;
  case 1: ret[0] = "stop on epsg"; break;
  case 2: ret[0] = "stop on epsf"; break;
  case 3: ret[0] = "stop on epsx"; break;
  case 4: ret[0] = "stop because of maximum function evaulations"; break;
  case 5: ret[0] = "stop because of maximum iterations"; break;
  case 6: ret[0] = "slope in the opposite direction to the gradient too small"; break;
  case 7: ret[0] = "stop when calculating the descent direction"; break;
  case 8: ret[0] = "stop when calculating the Hessian approximation"; break;
  case 10: ret[0] = "stop by failure of the linear search, cause not specified"; break;
  case 11: ret[0] = "stop by failure of the linear search, cause not specified with indsim <0"; break;
  case 12: ret[0] = "a step too small close to a step too big this can result from an error in the gradient"; break;
  case 13: ret[0] = "too many calls in a linear search"; break;
  default:
    ret[0] = "";
  }
  
  delete[] x;
  delete[] g;
  delete[] lower;
  delete[] upper;
  delete[] itrav;
  delete[] trav;
  return Rcpp::List::create(Rcpp::Named("value") = f,
			    Rcpp::Named("par") = par,
			    // Rcpp::Named("L") = L,
			    // Rcpp::Named("D") = D,
			    Rcpp::Named("trav") = travS,
			    // Rcpp::Named("zm")=zms,
			    // Rcpp::Named("c.hess") = hess,
			    Rcpp::Named("code") = indqn[0],
                            Rcpp::Named("codeDesc") = ret,
			    Rcpp::Named("n.fn") = n1qn1_calls,
			    Rcpp::Named("n.gr") = n1qn1_grads);
  END_RCPP
    }
