#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

typedef void (*S2_fp) (int *, int *, double *, double *, double *, int *, float *, double *);
extern void n1qn1c_(S2_fp simul, int n[], double x[], double f[], double g[], double var[], double eps[],
		   int mode[], int niter[], int nsim[], int imp[], double zm[], int izs[], 
		   float rzs[], double dzs[]);
void n1qn1cF(S2_fp simul, int n[], double x[], double f[], double g[], double var[], double eps[],
		   int mode[], int niter[], int nsim[], int imp[], int lp[], double zm[], int izs[], 
		   float rzs[], double dzs[]) {
  n1qn1c_(simul, n, x, f, g, var, eps, mode, niter, nsim, imp, zm, izs, rzs, dzs);
}

void n1qn1cF2(S2_fp simul, int n[], double x[], double f[], double g[], double var[], double eps[],
		   int mode[], int niter[], int nsim[], int imp[], double zm[], int izs[], 
		   float rzs[], double dzs[]) {
  n1qn1c_(simul, n, x, f, g, var, eps, mode, niter, nsim, imp, zm, izs, rzs, dzs);
}
/* .C calls */
extern SEXP n1qn1c_wrap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

void R_init_n1qn1c(DllInfo *dll)
{
  R_CallMethodDef callMethods[]  = {
    {"n1qn1c_wrap", (DL_FUNC) &n1qn1c_wrap, 13},
    {NULL, NULL, 0}
  };
  R_RegisterCCallable("n1qn1c","n1qn1cF", (DL_FUNC) n1qn1cF);
  R_RegisterCCallable("n1qn1c","n1qn1cF2", (DL_FUNC) n1qn1cF2);
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
} 
