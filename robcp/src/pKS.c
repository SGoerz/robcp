#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>

// Cumulative distribution function of the limit of the two-sided KS-statistic
/* References:
 * van Mulbregt, P. (2018). 
 * Computing the Cumulative Distribution Function and Quantiles of the One-sided 
 * Kolmogorov-Smirnov Statistic. arXiv preprint arXiv:1802.06966.
 */
static void KSdist(int n, double *x, double tol)
{
  double x_new, x_old, s, z, t;
  int i, k, k_max;
  
  k_max = (int) sqrt(2 - log(tol));
  
  for(i = 0; i < n; i++) 
  {
    if(x[i] <= 0)
    {
      x[i] = 0;
    } else if(x[i] < 1) 
    {    
      t = exp(- M_PI * M_PI / (8 * x[i] * x[i]));
      z = 1 + 8 * pow(t, k_max);
      s = M_SQRT_PI * M_SQRT2 / x[i];
      
      for(k = (k_max - 1); k >= 1; k--)
      {
        z = (1 + pow(t, 8 * k) * z);
      }
      
      x[i] = s * t * z;
    } else 
    {
      z = -2 * x[i] * x[i];
      s = -1;
      k = 1;
      x_old = 0;
      x_new = 0.5;
      tol /= 2;
      while(fabs(x_old - x_new) > tol) 
      {
        x_old = x_new;
        x_new += s * exp(z * k * k);
        s *= -1;
        k++;
      }
      x[i] = 2 * x_new;
    }
  }
}

// wrapper for KSdist
SEXP pKSdist(SEXP statistic, SEXP stol)
{
  int n = LENGTH(statistic);
  double tol = asReal(stol);
  SEXP ans = duplicate(statistic);
  KSdist(n, REAL(ans), tol);
  return ans;
}
