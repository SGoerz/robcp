#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>

//** Wilcoxon-Mann-Whitney test **//

double hFun(double x1, double x2)
{
  if(x1 < x2) return 0.5;
  else if(x1 > x2) return -0.5;
  return 0;
}

SEXP wilcox(SEXP X)
{
  int n = length(X);
  double *x = REAL(X);
  
  SEXP RES;
  PROTECT(RES = allocVector(REALSXP, n-1));
  double *res = REAL(RES);

  int i, j, k;
  double sum = 0;
  
  for(j = 1; j < n; j++)
  {
    sum += hFun(x[0], x[j]);
  }
  res[0] = fabs(sum) / pow(sqrt(n), 3);

  for(k = 1; k < n-1; k++)
  {
    for(i = 0; i < k; i++)
    {
      sum -= hFun(x[i], x[k]);
    }
    
    for(i = k+1; i < n; i++)
    {
      sum += hFun(x[k], x[i]);
    }
    
    res[k] = fabs(sum) / pow(sqrt(n), 3);
  }
  
  /*
  for(k = 1; k < n; k++)
  {
    sum = 0;
    for(i = 0; i <= k; i++)
    {
      for(j = k+1; j < n; j++)
      {
        sum += hFun(x[i], x[j]);
      }
    }
    if(fabs(sum) > max[0])
    {
      max[0] = fabs(sum);
      max[1] = k + 1;
    }
  }
  */
  
  UNPROTECT(1);
  return RES;
}