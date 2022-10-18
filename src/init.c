#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "n1qn1.h"

extern void n1qn1_(S2_fp simul, int n[], double x[], double f[], double g[], double var[], double eps[],
		   int mode[], int niter[], int nsim[], int imp[], double zm[], int izs[], 
		   float rzs[], double dzs[], int id[]);
void n1qn1F(S2_fp simul, int n[], double x[], double f[], double g[], double var[], double eps[],
		   int mode[], int niter[], int nsim[], int imp[], int lp[], double zm[], int izs[], 
		   float rzs[], double dzs[]) {
  int id = 0;
  n1qn1_(simul, n, x, f, g, var, eps, mode, niter, nsim, imp, zm, izs, rzs, dzs, &id);
}

void n1qn1F2(S2_fp simul, int n[], double x[], double f[], double g[], double var[], double eps[],
		   int mode[], int niter[], int nsim[], int imp[], double zm[], int izs[], 
		   float rzs[], double dzs[]) {
  int id = 0;
  n1qn1_(simul, n, x, f, g, var, eps, mode, niter, nsim, imp, zm, izs, rzs, dzs, &id);
}
/* .C calls */
extern SEXP n1qn1_wrap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

void R_init_n1qn1(DllInfo *dll)
{
  R_CallMethodDef callMethods[]  = {
    {"n1qn1_wrap", (DL_FUNC) &n1qn1_wrap, 13},
    {NULL, NULL, 0}
  };
  R_RegisterCCallable("n1qn1","n1qn1F", (DL_FUNC) &n1qn1F);
  R_RegisterCCallable("n1qn1","n1qn1F2", (DL_FUNC) &n1qn1F2);
  R_RegisterCCallable("n1qn1","n1qn1_", (DL_FUNC) &n1qn1_);
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
} 
