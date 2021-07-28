#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>

//** Wilcoxon-Mann-Whitney test **//

double h1(double x1, double x2)
{
  if(x1 < x2) return 0.5;
  else return -0.5;
}

double h2(double x1, double x2)
{
  return x1 - x2;
}

SEXP wilcox(SEXP X, SEXP H)
{
  int n = length(X);
  double *x = REAL(X);
  int h = *REAL(H);
  
  double (*hFun)(double, double);
  switch(h)
  {
  case 1: hFun = &h1; break;
  case 2: hFun = &h2; break;
  default: hFun = &h1; break;
  }
  
  SEXP MAX;
  PROTECT(MAX = allocVector(REALSXP, 2));
  double *max = REAL(MAX);
  max[0] = 0;

  int i, j, k;
  double sum = 0;
  
  for(j = 1; j < n; j++)
  {
    sum += hFun(x[0], x[j]);
  }
  max[0] = fabs(sum);
  max[1] = 1;
  
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
    
    if(fabs(sum) > max[0])
    {
      max[0] = fabs(sum);
      max[1] = k + 1;
    }
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
  
  max[0] = max[0] / pow(sqrt(n), 3);
  
  UNPROTECT(1);
  return MAX;
}