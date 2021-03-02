#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>

//** code for the computation of cumulative sums **//

/* h_cumsum: cumulative sum
 * 
 * input: Y (numeric vector)
 * 
 * output: cumulative sum (numeric vector)
 */
SEXP h_cumsum(SEXP Y)
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

/* cumsum_ma: columnwise cumulative sum of a matrix
 * 
 * input: Y (numeric 'matrix' in shape of a vector (along columns))
 *        N,M (number of rows and columns of Y)
 * 
 * output: numeric 'matrix' of columnwise cumulative sums
 */
SEXP h_cumsum_ma(SEXP Y, SEXP N, SEXP M)
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


//** computes the test statistic for the Huberized change point test **//

/* h_teststat: test statistic for a single time series
 * 
 * input: Y (time series; numeric vector)
 *        SIGMA (inverted estimated long run variance; numeric)
 * 
 * output: test statistic (numeric)
 */
SEXP h_teststat(SEXP Y, SEXP SIGMA)
{
  PROTECT(Y);
  double sigma = *REAL(SIGMA);
  
  SEXP MAX;
  PROTECT(MAX = allocVector(REALSXP, 1));
  double *max = REAL(MAX);
  max[0] = 0;
  
  int n = length(Y);
  
  double sqn = sqrt(n);
  // does that work?????????????????????????????????????????????????????????????
  double *cumsum = REAL(h_cumsum(Y));
  double sumN = cumsum[n - 1] / n;
  double temp;
  
  int i;
  
  for(i = 0; i < n; i++)
  {
    temp = fabs(cumsum[i] - (i + 1) * sumN);
    
    if(temp > max[0]) max[0] = temp; 
  }
  max[0] /= (sqn * sigma);
  
  UNPROTECT(2);
  return(MAX);
}


/* h_teststat_ma: test statistic for a more-dim. time series
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
 * output: test statistic (numeric) 
 */
SEXP h_teststat_ma(SEXP Y, SEXP SIGMA, SEXP SWAPS, SEXP N, SEXP M)
{
  PROTECT(Y);
  PROTECT(SIGMA);
  PROTECT(SWAPS);
  double *sigma = REAL(SIGMA);
  double *swaps = REAL(SWAPS);
  
  int n = *REAL(N);
  int m = *REAL(M);
  
  /* cumsum: matrix; each column is the cumulative sum vector of the 
   *         corresponding column in Y
   */
  // does that work?????????????????????????????????????????????????????????????
  double *cumsum = REAL(h_cumsum_ma(Y, N, M));
  
  double temp[m];
  double maxCand, temp2;
  
  SEXP MAX; 
  PROTECT(MAX = allocVector(REALSXP, 1));
  double *max = REAL(MAX);
  max[0] = 0;
  
  int i, j, k, index;
  
  for(i = 0; i < n; i++)
  {
    for(j = 0; j < m; j++)
    {
      // temp: CUSUM statistic for each column
      temp[j] = cumsum[j * n + i] - (i + 1) * cumsum[(j + 1) * n - 1] / n;
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
    
    maxCand = 0;
    
    for(j = 0; j < m; j++)
    {
      for(k = j; k < m; k++)
      {
        if(j == k)
        {
          maxCand += temp[j]*temp[j] * sigma[j + j * m];
        }
        else
        {
          maxCand += 2 * (temp[j]*temp[k] * sigma[k + j * m]);
        }
      }
    }
    maxCand /= n;
    
    if(maxCand > max[0]) max[0] = maxCand;
  }
  
  UNPROTECT(4);
  return MAX;
}