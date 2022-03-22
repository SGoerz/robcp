#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>

int comp(const void *elem1, const void *elem2) 
{
  int f = *((int*)elem1);
  int s = *((int*)elem2);
  if (f > s) return  1;
  if (f < s) return -1;
  return 0;
}

//** code for the computation of cumulative sums **//

/* c_cumsum: cumulative sum
 * 
 * input: Y (numeric vector)
 * 
 * output: cumulative sum (numeric vector)
 */
SEXP c_cumsum(SEXP Y)
{
  SEXP X = duplicate(Y);
  PROTECT(X);
  double *x = REAL(X);
  
  int n = length(X);
  int i;
  
  for(i = 1; i < n; i++)
  {
    x[i] += x[i - 1];
  }
  UNPROTECT(1);
  return(X);
}

/* c_cumsum_ma: columnwise cumulative sum of a matrix
 * 
 * input: Y (numeric 'matrix' in shape of a vector (along columns))
 *        N,M (number of rows and columns of Y)
 * 
 * output: numeric 'matrix' of columnwise cumulative sums
 */
SEXP c_cumsum_ma(SEXP Y, SEXP N, SEXP M)
{
  int n = *REAL(N);
  int m = *REAL(M);
  
  SEXP X = duplicate(Y);
  PROTECT(X);
  double *x = REAL(X);
  
  int i, j;
  
  for(j = 0; j < m; j++)
  {
    for(i = 1; i < n; i++)
    {
      x[i + n * j] += x[i + n * j - 1];
    }
  }
  
  UNPROTECT(1);
  return X;
}

//** computes the test statistic for the CUSUM change point test **//

/* CUSUM: test statistic for a single time series
 * 
 * input: Y (time series; numeric vector)
 * 
 * output: test statistic (numeric vector)
 */
SEXP CUSUM(SEXP Y)
{
  PROTECT(Y);
  int n = length(Y);

  SEXP RES;
  PROTECT(RES = allocVector(REALSXP, n-1));
  double *res = REAL(RES);
  double sqn = sqrt(n);
  double *csum = REAL(c_cumsum(Y));
  double sumN = csum[n - 1] / n;
  
  int i;
  
  for(i = 0; i < n-1; i++)
  {
    res[i] = fabs(csum[i] - (i + 1) * sumN) / sqn;
  }

  UNPROTECT(2);
  return(RES);
}


/* CUSUM_ma: test statistic for a more-dim. time series
 * 
 * input: Y (time series; numeric 'matrix' in the form of a vector; columnwise)
 *        SIGMA (inverted estimated long run covariance; numeric 'matrix' in 
 *               form of a vector; columnwise)
 *        SWAPS (indices of the rows and columns that were swapped in x in order 
 *               to compute the modified Cholesky factorization. For example if 
 *               the i-th element of swaps is the number j, then the i-th and 
 *               the j-th row and column were swapped. To reconstruct the
 *               original matrix swaps has to be read backwards.)
 *        N (length of time series; integer)
 *        M (dimension of time series; integer)
 *        
 * output: test statistic (numeric vector) 
 */
SEXP CUSUM_ma(SEXP Y, SEXP SIGMA, SEXP SWAPS, SEXP N, SEXP M)
{
  PROTECT(Y);
  PROTECT(SIGMA);
  PROTECT(SWAPS);
  double *sigma = REAL(SIGMA);
  double *swaps = REAL(SWAPS);
  
  int n = *REAL(N);
  int m = *REAL(M);
  
  /* csum: matrix; each column is the cumulative sum vector of the 
   *       corresponding column in Y
   */
  double *csum = REAL(c_cumsum_ma(Y, N, M));
  
  double temp[m];
  double temp2;
  
  SEXP RES; 
  PROTECT(RES = allocVector(REALSXP, n));
  double *res = REAL(RES);
  
  int i, j, k, index;
  
  for(i = 0; i < n; i++)
  {
    for(j = 0; j < m; j++)
    {
      // temp: CUSUM statistic for each column
      temp[j] = csum[j * n + i] - (i + 1) * csum[(j + 1) * n - 1] / n;
    }

    // swap temp according to swaps
    for(j = 0; j < m; j++)
    {
      if(j != swaps[j])
      {
        index = swaps[j];
        
        temp2 = temp[j];
        temp[j] = temp[index];
        temp[index] = temp2;
      }
    }
    
    res[i] = 0;

    for(j = 0; j < m; j++)
    {
      for(k = j; k < m; k++)
      {
        if(j == k)
        {
          res[i] += temp[j]*temp[j] * sigma[j + j * m];
        }
        else
        {
          res[i] += 2 * (temp[j]*temp[k] * sigma[k + j * m]);
        }
      }
    }
    res[i] /= n;
  }
  
  UNPROTECT(4);
  return RES;
}


SEXP MD(SEXP X, SEXP CUMMED, SEXP N)
{
  double n = *REAL(N);
  double *x = REAL(X);
  double *cummed = REAL(CUMMED);
  
  SEXP RES; 
  PROTECT(RES = allocVector(REALSXP, n-1));
  double *res = REAL(RES);

  int i, k;
  
  for(k = 1; k < n; k++)
  {
    res[k-1] = 0;
    for(i = 0; i <= k; i++)
    {
      res[k-1] += fabs(x[i] - cummed[k]);
    }
  }
  
  UNPROTECT(1);
  return RES;
}


SEXP GMD(SEXP X, SEXP N)
{
  double n = *REAL(N);
  double *x = REAL(X);

  SEXP RES; 
  PROTECT(RES = allocVector(REALSXP, n-1));
  double *res = REAL(RES);
  
  int i, k;
  
  res[0] = fabs(x[0] - x[1]);
  
  for(k = 2; k < n; k++)
  {
    res[k-1] = res[k-2];
    for(i = 0; i < k; i++)
    {
      res[k-1] += fabs(x[i] - x[k]);
    }
  }
  
  UNPROTECT(1);
  return RES;
}