#include <R.h>
#include <Rinternals.h>

// modified from scilab to use R interface.  
int F77_SUB(vfinite)(int *n, double *v)
{
  for (int i = 0; i < *n; ++i) {
    if (!R_FINITE(v[i])) {
      return 0;
    }
  }
  return 1;
}


int F77_SUB(vfinite1)(double *v)
{
  if (!R_FINITE(v[0])) {
      return 0;
  }
  return 1;
}
