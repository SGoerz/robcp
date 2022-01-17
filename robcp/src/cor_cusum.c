#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>


SEXP tau(SEXP X, SEXP Y, SEXP N)
{
  double *x = REAL(X);
  int n = *REAL(N);
  double *y = REAL(Y);

  SEXP ERG;
  PROTECT(ERG = allocVector(REALSXP, n-1));
  double *erg = REAL(ERG);
  
  int i, j;
  double prd, count;
  
  count = 0;
  
  for(j = 1; j < n; j++)
  {
    for(i = 0; i < j; i++)
    {
      prd = (x[j] - x[i]) * (y[j] - y[i]);
      if(prd < 0) count -= 1; else if(prd > 0) count += 1;
    }
    erg[j-1] = count * 2 / (j * (j+1));
  }
  
  UNPROTECT(1);
  return ERG;
}