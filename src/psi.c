#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>

//** code for standardizing and psi-transforming the time series **//


// standardizes arr[start:(n-1)] by location mu and scale sigma
void trafo(double arr[], double mu, double sigma, int start, int n)
{
  int i;
  
  for(i = start; i < start + n; i++)
  {
    arr[i] = (arr[i] - mu) / sigma;
  }
}


// marginal Huberized location
void HLm(double arr[], int start, int n, int m, double k)
{
  int j; 
  for(j = 0; j < m; j++)
  {
    if(arr[start + j * n] > k) arr[start + j * n] = k;
    else if(arr[start + j * n] < -k) arr[start + j * n] = -k;
  }
}

// global Huberized location
void HLg(double arr[], int start, int n, int m, double k)
{
  int j; 
  double sum = 0;
  
  for(j = 0; j < m; j++)
  {
    sum = sum + (arr[start + j * n] * arr[start + j * n]);
  }
  sum = sqrt(sum);
  
  if(sum == 0)
  {
    for(j = 0; j < m; j++)
    {
      arr[start + j * n] = 0;
    }
  }
  else if(sum > k)
  {
    for(j = 0; j < m; j++)
    {
      arr[start + j * n] = arr[start + j * n] * k / sum;
    }
  }
}


// marginal sign location
void SLm(double arr[], int start, int n, int m)
{
  int j;
  for(j = 0; j < m; j++)
  {
    if(arr[start + j * n] < 0) arr[start + j * n] = -1;
    else if(arr[start + j * n] > 0) arr[start + j * n] = 1;
  }
}


// global sign location
void SLg(double arr[], int start, int n, int m, double k)
{
  int j; 
  double sum = 0;
  
  for(j = 0; j < m; j++)
  {
    sum = sum + (arr[start + j * n] * arr[start + j * n]);
  }
  sum = sqrt(sum);
  
  if(sum == 0)
  {
    for(j = 0; j < m; j++)
    {
      arr[start + j * n] = 0;
    }
  }
  else 
  {
    for(j = 0; j < m; j++)
    {
      arr[start + j * n] = arr[start + j * n] / sum;
    }
  }
}


/* psi_location: standartizes and transforms a time series Y according to a psi 
 *               function FUN for the location
 *               
 * input: Y (time series; numeric 'matrix')
 *        FUN (psi function; integer)
 *        N (nrow(Y); integer)
 *        M (ncol(Y); integer)
 *        K (parameter for psi function; numeric)
 *        MED (columnwise median of Y; numeric vector)
 *        MAD (columnwise mad of Y; numeric vector)
 *        
 * output: transformed time series (numeric 'matrix')
 */
SEXP psi_location(SEXP Y, SEXP FUN, SEXP N, SEXP M, SEXP K, //SEXP CONST, 
                  SEXP MED, SEXP MAD)
{
  int n = *REAL(N);
  int m = *REAL(M);
  int fun = *REAL(FUN);
  //double c = *REAL(CONST);
  double k = *REAL(K);
  double *med = REAL(MED);
  double *mad = REAL(MAD);
  // double med, mad;
  
  SEXP X = duplicate(Y);
  PROTECT(X);
  double *x = REAL(X);
  
  int i, j;
  //double temp[n];
  
  for(j = 0; j < m; j++)
  {
    //// copy columns of x (matrix) into temp 
    // extract(x, temp, j * n, n);
    
    //// for now: median and mad need to be computed in R
    //med = median(temp, n);
    //mad = fastmad(temp, med, c, n);
    
    trafo(x, med[j], mad[j], j * n, n);
  }
  for(i = 0; i < n; i++)
  {
    // apply psi function
    switch(fun)
    {
      case 1: HLm(x, i, n, m, k); break;
      case 2: HLg(x, i, n, m, k); break;
      case 3: SLm(x, i, n, m);    break;
      case 4: SLg(x, i, n, m, k); break;
    }
  }
  
  UNPROTECT(1);
  return X;
}

/* psi_covariance: standartizes and transforms a time series Y according to a 
 *                 psi function FUN for the covariance
 *                 
 * input: Y (time series; numeric 'matrix')
 *        FUN (psi function; integer)
 *        N (nrow(Y); integer)
 *        M (ncol(Y); integer)
 *        K (parameter for psi function; numeric)
 *        MED (columnwise median of Y; numeric vector)
 *        MAD (columnwise mad of Y; numeric vector)
 *        
 * output: transformed time series (numeric 'matrix')
 */
SEXP psi_covariance(SEXP Y, SEXP FUN, SEXP N, SEXP M, SEXP K, //SEXP CONST, 
                    SEXP MED, SEXP MAD)
{
  int fun = *REAL(FUN);
  SEXP FUN2;
  PROTECT(FUN2 = allocVector(REALSXP, 1));
  double *fun2 = REAL(FUN2);
  fun2[0] = fun - 4;
  // FUN2 now indicates the corresponding location transformation
  
  // apply location tranformation
  SEXP X = psi_location(Y, FUN2, N, M, K, //CONST,
                        MED, MAD);
  PROTECT(X);
  double *x = REAL(X);
  
  int n = *REAL(N);
  int m = *REAL(M);
  
  // dimension of the result
  int dimErg;     
  switch(fun)
  {
    case 7: dimErg = m * (m - 1) / 2;     break;
    case 8: dimErg = m * (m + 1) / 2 - 1; break; 
    default: dimErg = m * (m + 1) / 2;    break;
  }    
  SEXP ERG;
  PROTECT(ERG = allocVector(REALSXP, n * dimErg));
  double *erg = REAL(ERG);
  
  int i, j, k, l, upper;
  
  if(fun == 8) upper = m - 1;
  else upper = m;
  
  for(i = 0; i < n; i++)
  {
    l = 0;
    
    if(fun == 7)
    {
      for(j = 0; j < upper; j++)
      {
        for(k = (j + 1); k < m; k++)
        {
          erg[i + l * n] = x[i + j * n] * x[i + k * n];
          l++;
        }
      }
    }
    else
    {
      for(j = 0; j < upper; j++)
      {
        for(k = j; k < m; k++)
        {
          erg[i + l * n] = x[i + j * n] * x[i + k * n];
          l++;
        }
      }
    }
  }
  
  UNPROTECT(3);
  return ERG;
}


// psi: wrapper for psi_location and psi_covariance
SEXP psi(SEXP Y, SEXP FUN, SEXP N, SEXP M, SEXP K, //SEXP CONST, 
         SEXP MED, SEXP MAD)
{
  int fun = *REAL(FUN);
  
  SEXP ERG;
  if(fun <= 4) ERG = psi_location(Y, FUN, N, M, K, //CONST, 
     MED, MAD);
  else ERG = psi_covariance(Y, FUN, N, M, K, //CONST, 
                            MED, MAD);
  
  return ERG;
}