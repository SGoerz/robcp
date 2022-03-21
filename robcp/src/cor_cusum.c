#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>

/* tau: estimate Kendall's tau for x[1:k], k = 2, ..., n
 * 
 * input: X (first part of the bivariate time series (numeric vector))
 *        X (second part of the bivariate time series (numeric vector))
 *        n (length of the time series; numeric or integer)
 *        
 * output: vector of tau-s (numeric; length n-1)
 */
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