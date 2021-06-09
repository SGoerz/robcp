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


//*************************** kernel functions *******************************//

// Bartlett kernel 
double kBartlett(double x)
{
  if(fabs(x) < 1) return 1 - fabs(x);
  else return 0;
}

// flat top kernel  
double kFT(double x)
{
  if(fabs(x) <= 0.5) return 1;
  else if(fabs(x) > 0.5 && fabs(x) < 1) return (2 - 2 * fabs(x));
  else return 0;
}

// Parzen kernel
double kParzen(double x)
{
  if(0 <= fabs(x) && fabs(x) <= 0.5)
  {
    return 1 - 6 * x * x + 6 * x * x * fabs(x);
  } else if(0.5 < fabs(x) && fabs(x) <= 1)
  {
    return 2 * pow(1 - fabs(x), 3);
  }
  return 0;
}

// Quadratic Spectral kernel
double kQS(double x)
{
  if(x == 0) return 1;
  return 25 / (12 * M_PI * M_PI * x * x) *
    (sin(6 * M_PI * x / 5) / (6 * M_PI* x / 5) - cos(6 * M_PI * x / 5)); 
}

// Tukey-Hanning kernel
double kTH(double x)
{
  if(fabs(x) <= 1) return (1 + cos(M_PI * x)) / 2;
  return 0;
}

// truncated kernel
double kTruncated(double x)
{
  if(fabs(x) > 1) return 0;
  return 1;
}


//*************** kernel-bases long run variance estimation ******************//

/* sigma_1: estimate the long run variance of a single signal
 * 
 * input: x (numeric centered vector)
 *        n (length of x)
 *        b_n (bandwidth; numeric or integer)
 *        
 * output: estimated long run variance (numeric)
 */
double sigma_1(double *x, int n, double b_n, int k)
{
  int i, h;
  
  double (*kFun)(double);
  
  switch(k)
  {
  case 1: kFun = &kBartlett; break;
  case 2: kFun = &kFT; break;
  case 3: kFun = &kParzen; break;
  case 4: kFun = &kQS; break;
  case 5: kFun = &kTH; break;
  case 6: kFun = &kTruncated; break;
  default: kFun = &kBartlett; break;
  }
  
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
    temp += temp2 * kFun(h / b_n);
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
double sigma_2(double x1[], double x2[], int n, double b_n, int k)
{
  int i, h;
  
  double (*kFun)(double);
  
  switch(k)
  {
  case 1: kFun = &kBartlett; break;
  case 2: kFun = &kFT; break;
  case 3: kFun = &kParzen; break;
  case 4: kFun = &kQS; break;
  case 5: kFun = &kTH; break;
  case 6: kFun = &kTruncated; break;
  default: kFun = &kBartlett; break;
  }
  
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
    erg += temp * kFun(h / b_n);
  }
  
  erg /= n;
  
  return erg;
}

// wrapper for the long run variance function sigma_1 in the one-dim. case
// long run variance estimation via kernel density
SEXP lrv(SEXP X, SEXP BN, SEXP K)
{
  double *x = REAL(X);
  int n = length(X);
  double b_n = *REAL(BN);
  int k = *REAL(K);
  
  SEXP ERG;
  PROTECT(ERG = allocVector(REALSXP, 1));
  
  double *erg = REAL(ERG);
  erg[0] = sigma_1(x, n, b_n, k);
  
  UNPROTECT(1);
  return ERG;
}

/* wrapper for the long run covariance matrix using sigma_1 (variance) and 
 * sigma_2 (covariance)
 */
SEXP lrv_matrix(SEXP Y, SEXP N, SEXP M, SEXP BN, SEXP K)
{
  SEXP X = duplicate(Y);
  PROTECT(X);
  double *x = REAL(X);
  
  int n = *REAL(N);
  int m = *REAL(M);
  double b_n = *REAL(BN);
  int k = *REAL(K);
  
  SEXP ERG;
  PROTECT(ERG = allocVector(REALSXP, m * m));
  double *erg = REAL(ERG);
  
  int j, i;
  
  double arr1[n];
  double arr2[n];
  double temp;
  
  for(j = 0; j < m; j++)
  {
    for(i = j; i < m; i++)
    {
      if(j == i)
      {                      
        extract(x, arr1, j * n, n);
        erg[j + j * m] = sigma_1(arr1, n, b_n, k);
      }
      else
      {
        extract(x, arr1, j * n, n);
        extract(x, arr2, i * n, n);
        
        temp = sigma_2(arr1, arr2, n, b_n, k);
        
        erg[j + i * m] = temp;
        erg[i + j * m] = temp;
      }
    }
  }
  
  UNPROTECT(2);
  return ERG;
}


//********************* for dependent wild bootstrap *************************//

// gen_matrix: generate matrix for the dependent wild bootstrap
//
// input: N (sample size)
//        L (block size)
//
// output: NxN 'matrix'; columnwise as a vector
SEXP gen_matrix(SEXP N, SEXP L, SEXP K)
{
  int n = *REAL(N);
  int l = *REAL(L);
  int k = *REAL(K); 
  double (*kFun)(double);
  
  switch(k)
  {
  case 1: kFun = &kBartlett; break;
  case 3: kFun = &kParzen; break;
  case 4: kFun = &kQS; break;
  default: kFun = &kBartlett; break;
  }
  
  SEXP MATRIX; 
  PROTECT(MATRIX = allocVector(REALSXP, n*n));
  double *matrix = REAL(MATRIX);
  
  int i, j;
  
  for(i = 0; i < n; i++)
  {
    for(j = i; j < n; j++)
    {
      matrix[i * n + j] = kFun((double)(i - j) / l);
      if(i != j) matrix[j * n + i] = matrix[i * n + j];
    }
  }
  
  UNPROTECT(1);
  return MATRIX;
}


//********************** other variance estimation ***************************//


//// u
SEXP u_hat(SEXP X, SEXP B, SEXP K)
{
  SEXP SUM;
  PROTECT(SUM = allocVector(REALSXP, 1));
  double *sum = REAL(SUM);
  sum[0] = 0;
  
  double *x = REAL(X);
  double b = *REAL(B);
  int n = length(X);
  int k = *REAL(K);
  double (*kFun)(double);
  
  switch(k)
  {
  case 1: kFun = &kBartlett; break;
  case 2: kFun = &kFT; break;
  case 3: kFun = &kParzen; break;
  case 4: kFun = &kQS; break;
  case 5: kFun = &kTH; break;
  case 6: kFun = &kTruncated; break;
  default: kFun = &kBartlett; break;
  }
  
  int i, j;
  
  for(i = 0; i < n; i++)
  {
    for(j = i+1; j < n; j++)
    {
      sum[0] += kFun((x[i] - x[j]) / b);
    }
  }
  
  sum[0] = sum[0] * 2 / (n * (n - 1) * b);
  
  UNPROTECT(1);
  return SUM;
}

// lrv_subs_nonoverlap: non-overlapping subsampling estimation of the long run variance
//
// input: X: empirical cumulative distribution function
//        L: block length
//        MEAN: mean of X * L
//        DISTR: plain observations or emp. distribution function?
SEXP lrv_subs_nonoverlap(SEXP X, SEXP L, SEXP MEAN, SEXP DISTR)
{
  SEXP SUM;
  PROTECT(SUM = allocVector(REALSXP, 1));
  double *sum = REAL(SUM);
  sum[0] = 0;
  
  double *x = REAL(X);
  int l = *REAL(L);
  double mean = *REAL(MEAN);
  int distr = *REAL(DISTR);
  int n = length(X);
  
  int i, j;
  double temp;
  
  double blocknr = n / l; // ??? wie Runden ???

  for(i = 1; i <= blocknr; i++)
  {
    temp = 0;
    
    for(j = (i - 1) * l; j < i * l; j++)
    {
      temp += x[j];
    }
    
    temp = temp - mean;
    if(distr == 1)
    {
      sum[0] += fabs(temp);
    } else
    {
      sum[0] += temp * temp;
    }
  }
  
  sum[0] /= blocknr;
  
  if(distr == 1)
  {
    sum[0] = sum[0] * sqrt(M_PI_2 / l);
  } else
  {
    sum[0] = sum[0] / l;
  }
  
  UNPROTECT(1);
  return SUM;
}

// lrv_subs_overlap: overlapping subsampling estimation of the long run variance
//
// input: X: empirical cumulative distribution function
//        L: block length
SEXP lrv_subs_overlap(SEXP X, SEXP L)
{
  SEXP SUM;
  PROTECT(SUM = allocVector(REALSXP, 1));
  double *sum = REAL(SUM);
  sum[0] = 0;
  
  double *x = REAL(X);
  int l = *REAL(L);
  int n = length(X);
  
  int i, j;
  double temp;
  
  temp = 0;
  for(j = 0; j < l; j++)
  {
    temp += x[j];
  }
  sum[0] += fabs(temp - l * 0.5);
  
  for(i = 1; i <= n - l; i++)
  {
    temp = temp - x[i - 1] + x[i + l - 1];
    sum[0] += fabs(temp - l * 0.5);
  }
  
  sum[0] = sum[0] * sqrt(M_PI) / (sqrt(2 * l) * (n - l + 1));

  UNPROTECT(1);
  return SUM;
}
