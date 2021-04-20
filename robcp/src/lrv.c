#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>

//** estimation of the long run variance resp. covariance matrix **//


//function to extract a column of a matrix (given by vector ma[]) into arr[]
void extract(double ma[], double arr[], int start, int n)
{
  int i;
  
  for(i = 0; i < n; i++)
  {
    arr[i] = ma[start + i];
  }
}


// flat top kernel function 
double k(double x)
{
  if(fabs(x) <= 0.5) return 1;
  else if(fabs(x) > 0.5 && fabs(x) < 1) return (2 - 2 * fabs(x));
  else return 0;
}


/* sigma_1: estimate the long run variance of a single signal
 * 
 * input: x (numeric centered vector)
 *        n (length of x)
 *        b_n (bandwidth; numeric or integer)
 *        
 * output: estimated long run variance (numeric)
 */
double sigma_1(double *x, int n, double b_n)
{
  int i, h;
  
  double var = 0;
  double temp = 0;
  double temp2;
  double erg = 0;
  
  for(i = 0; i < n; i++)
  {
    var += x[i] * x[i];
  }
  
  for(h = 1; h <= b_n; h++)
  {
    temp2 = 0;
    for(i = 0; i < (n - h); i++)
    {
      temp2 += x[i] * x[i + h];
    }
    temp += temp2 * k(h / b_n);
  }
  
  erg = (var + 2 * temp) / n;
  return erg;
}

/* sigma_2: long run covariance of two signals
 * 
 * input: x1, x2 (time series; numeric vectors)
 *        n (length of the time series; integer)
 *        b_n (bandwidth; numeric or integer)
 *        
 * output: long run covariance (numeric) 
 */
double sigma_2(double x1[], double x2[], int n, double b_n)
{
  int i, h;
  
  double erg = 0;
  double temp;
  
  for(i = 0; i < n; i++)
  {
    erg += (x1[i] * x2[i]);
  }
  for(h = 1; h <= b_n; h++)
  {
    temp = 0;
    for(i = 0; i < (n - h); i++)
    {
      temp += (x1[i] * x2[i + h] + x2[i] * x1[i + h]);
    }
    erg += temp * k(h / b_n);
  }
  
  erg /= n;
  
  return erg;
}

// wrapper for the long run variance function sigma_1 in the one-dim. case
SEXP lrv(SEXP X, SEXP BN)
{
  double *x = REAL(X);
  int n = length(X);
  double b_n = *REAL(BN);
  
  SEXP ERG;
  PROTECT(ERG = allocVector(REALSXP, 1));
  
  double *erg = REAL(ERG);
  erg[0] = sigma_1(x, n, b_n);
  
  UNPROTECT(1);
  return ERG;
}

/* wrapper for the long run covariance matrix using sigma_1 (variance) and 
 * sigma_2 (covariance)
 */
SEXP lrv_matrix(SEXP Y, SEXP N, SEXP M, SEXP BN)
{
  SEXP X = duplicate(Y);
  PROTECT(X);
  double *x = REAL(X);
  
  int n = *REAL(N);
  int m = *REAL(M);
  double b_n = *REAL(BN);
  
  SEXP ERG;
  PROTECT(ERG = allocVector(REALSXP, m * m));
  double *erg = REAL(ERG);
  
  int j, k;
  
  double arr1[n];
  double arr2[n];
  double temp;
  
  for(j = 0; j < m; j++)
  {
    for(k = j; k < m; k++)
    {
      if(j == k)
      {                      
        extract(x, arr1, j * n, n);
        erg[j + j * m] = sigma_1(arr1, n, b_n);
      }
      else
      {
        extract(x, arr1, j * n, n);
        extract(x, arr2, k * n, n);
        
        temp = sigma_2(arr1, arr2, n, b_n);
        
        erg[j + k * m] = temp;
        erg[k + j * m] = temp;
      }
    }
  }
  
  UNPROTECT(2);
  return ERG;
}





//// u
SEXP u_hat(SEXP X, SEXP B)
{
  SEXP SUM;
  PROTECT(SUM = allocVector(REALSXP, 1));
  double *sum = REAL(SUM);
  sum[0] = 0;
  
  double *x = REAL(X);
  double b = *REAL(B);
  int n = length(X);
  
  int i, j;
  
  for(i = 0; i < n; i++)
  {
    for(j = i+1; j < n; j++)
    {
      sum[0] += k((x[i] - x[j]) / b);
    }
  }
  
  sum[0] = sum[0] * 2 / (n * (n - 1) * b);
  
  UNPROTECT(1);
  return SUM;
}

///// sigma_hat_HL
SEXP lrv_subs(SEXP ECDF, SEXP L)
{
  SEXP SUM;
  PROTECT(SUM = allocVector(REALSXP, 1));
  double *sum = REAL(SUM);
  sum[0] = 0;
  
  double *ecdf = REAL(ECDF);
  double l = *REAL(L);
  int n = length(ECDF);
  
  int i, j;
  double temp;
  for(i = 0; i <= n - l; i++)
  {
    temp = 0;
    for(j = i; j < i + l; j++)
    {
      temp += ecdf[j];
    }
    sum[0] += fabs(temp - l * 0.5);
  }
  
  sum[0] = sum[0] * sqrt(M_PI) / (sqrt(2 * l) * (n - l + 1));
  
  UNPROTECT(1);
  return SUM;
}


SEXP gen_matrix(SEXP N, SEXP L)
{
  int n = *REAL(N);
  int l = *REAL(L);
  
  SEXP MATRIX; 
  PROTECT(MATRIX = allocVector(REALSXP, n*n));
  double *matrix = REAL(MATRIX);
  
  int i, j;
  
  for(i = 0; i < n; i++)
  {
    for(j = i; j < n; j++)
    {
      matrix[i * n + j] = k((double)(i - j) / l);
      if(i != j) matrix[j * n + i] = matrix[i * n + j];
    }
  }
  
  UNPROTECT(1);
  return MATRIX;
}
