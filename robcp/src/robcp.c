#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>

// k function for 
double k(double x)
{
       if(fabs(x) <= 0.5) return x;
       else if(fabs(x) > 0.5 && fabs(x) < 1) return (2 - 2 * fabs(x));
       else return 0;
}


//function to extract a column of a matrix (given by vector arr[]) into temp
void extract(double ma[], double arr[], int start, int n)
{
       int i;
       
       for(i = 0; i < n; i++)
       {
             arr[i] = ma[start + i];
       }
}


// sigma_1: long run variance of a single signal
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


// sigma_2: long run variance of two signals
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

// return sigma^2 for the one-dimensional case
SEXP sigma2(SEXP X, SEXP BN)
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


SEXP sigma_matrix(SEXP Y, SEXP N, SEXP M, SEXP BN)
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
                      erg[j + k * m] = sigma_1(arr1, n, b_n);
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

// finds the minimum in array arr with boundaries l and u
double minimum(double arr[], int l, int u)
{
       double min = arr[l];
       int i;
       
       for(i = l + 1; i <= u; i++)
       {
             if(arr[i] < min){min = arr[i];}
       }
       return min;
}

// standardizes arr[start:(n-1)] by location mu and scale sigma
void trafo(double arr[], double mu, double sigma, int start, int n)
{
     int i;
          
     for(i = start; i < start + n; i++)
     {
           arr[i] = (arr[i] - mu) / sigma;
     }
}


// marginal huberized location
void HLm(double arr[], int start, int n, int m, double k)
{
     int j; 
     for(j = 0; j < m; j++)
     {
           if(arr[start + j * n] > k) arr[start + j * n] = k;
           else if(arr[start + j * n] < -k) arr[start + j * n] = -k;
     }
}


// global huberized location
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
            for(j = 0; j < n; j++)
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
void VLm(double arr[], int start, int n, int m)
{
     int j;
     for(j = 0; j < m; j++)
     {
           if(arr[start + j * n] < 0) arr[start + j * n] = -1;
           else if(arr[start + j * n] > 0) arr[start + j * n] = 1;
     }
}


// global sign location
void VLg(double arr[], int start, int n, int m, double k)
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
            for(j = 0; j < n; j++)
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


SEXP psi_location(SEXP Y, SEXP FUN, SEXP N, SEXP M, SEXP K, //SEXP CONST, 
                  SEXP MED, SEXP MAD)
{
     int n = *REAL(N);
     int m = *REAL(M);
     int fun = *REAL(FUN);
     //double c = *REAL(CONST);
     double k = *REAL(K);
     
     SEXP X = duplicate(Y);
     PROTECT(X);
     double *x = REAL(X);
     
     int i, j;
     double temp[n];
     
     double *med = REAL(MED);
     double *mad = REAL(MAD);
     // double med, mad;
     
     for(j = 0; j < m; j++)
     {
           // copy parts of x (columns of matrix) into temp 
           extract(x, temp, j * n, n);
           
           // for now: median and mad need to be computed in R
           //med = median(temp, n);
           //mad = fastmad(temp, med, c, n);
           
           trafo(x, med[j], mad[j], j * n, n);
     }
     for(i = 0; i < n; i++)
     {
           switch(fun)
           {
                  case 1: HLm(x, i, n, m, k); break;
                  case 2: HLg(x, i, n, m, k); break;
                  case 3: VLm(x, i, n, m);    break;
                  case 4: VLg(x, i, n, m, k); break;
           }
     }
     
     UNPROTECT(1);
     return X;
}

SEXP psi_covariance(SEXP Y, SEXP FUN, SEXP N, SEXP M, SEXP K, //SEXP CONST, 
                    SEXP MED, SEXP MAD)
{
     int fun = *REAL(FUN);
     SEXP FUN2;
     PROTECT(FUN2 = allocVector(REALSXP, 1));
     double *fun2 = REAL(FUN2);
     fun2[0] = fun - 4;
     // FUN2 now indicates the corresponding location transformation
     
     SEXP X = psi_location(Y, FUN2, N, M, K, //CONST,
                           MED, MAD);
     PROTECT(X);
     double *x = REAL(X);
     
     int n = *REAL(N);
     int m = *REAL(M);
     
     int dimErg;     
     switch(fun)
     {
            case 7: dimErg = m * (m - 1) / 2; break;
            case 8: dimErg = m * (m + 1) / 2 - 1; break; 
            default: dimErg = m * (m + 1) / 2; break;
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


SEXP cumsum_ma(SEXP Y, SEXP N, SEXP M)
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


// Following: Revised Modified Cholesky Decomposition Algorithm

//swap rows and columns number j and maxindex
void rowColSwap(double A[], int j, int maxindex, int n)
{
     int i;
     double temp;
     
     for(i = 0; i < n; i++)
     {
           //Rows
           temp = A[i + j * n];
           A[i + j * n] = A[i + maxindex * n];
           A[i + maxindex * n] = temp;
     }
     
     for(i = 0; i < n; i++)
     {           
           //Columns
           temp = A[j + i * n];
           A[j + i * n] = A[maxindex + i * n];
           A[maxindex + i * n] = temp;
     }
}


void jthFac(double A[], double L[], int j, int n)
{
     int i, k;
     
     L[j + j * n] = sqrt(A[j + j * n]); 

     for(i = (j + 1); i < n; i++)
     {
           L[j + i * n] = A[j + i * n] / L[j + j * n]; 
           for(k = (j + 1); k < (i + 1); k++)       
           {
                 A[k + i * n] -= L[j + i * n] * L[j + k * n];
           }
     }
     for(i = (j + 1); i < n; i++)  

     {
           for(k = (i + 1); k < n; k++) 
           {
                A[k + i * n] = A[i + k * n];
           }
     }
}

//Revised Modified Cholesky Decomposition Algorithm
//A symmetric, stored in lower triangle
//tau = (macheps)^(1/3), tau.bar = (macheps)^(2/3), mu = 0.1
//goal: find LL^T of A + E, E >= 0
// Source of peusocode:
// Schnabel, R. B., & Eskow, E. (1999). "A revised modified Cholesky factorization algorithm" SIAM Journal on optimization, 9(4), 1135-1148.
void RMCDA(double A[], double L[], int n, double tau, double tau_bar, double mu, 
           double swaps[])
{
       int phaseone = 1; // boolean variable
       
       double gamma = fabs(A[0]);
       int i;
       int k = 0; 
       int maxindex = 0;
       double max = 0;
       double min = 0;
       double temp;
       for (i = 0; i < n * n - 1; i++) 
       { 
           L[i] = 0;
       }       
       //gamma = max(|A_ii|)
       for(i = 1; i < n; i++)
       {
             if(fabs(A[i + i * n]) > gamma) gamma = fabs(A[i + i * n]);
       }
       int j = 0;
       
       //Phase one, A potentially positive definite
       while(j < n && phaseone == 1)
       {
               max = A[j + j * n];
               maxindex = j;
               min = max;
               
               for(i = (j + 1); i < n; i++)
               {
                     if(A[i + i * n] > max) 
                     {
                            max = A[i + i * n];
                            maxindex = i;
                     }
                     if(A[i + i * n] < min)
                     {
                            min = A[i + i * n];
                     }
               }
               if((max < gamma * tau_bar) || (min < - mu * max))
               {
                       phaseone = 0;
               }
               else
               {
                       if(maxindex != j)
                       {
                                   rowColSwap(A, j, maxindex, n);
                                   rowColSwap(L, j, maxindex, n);
                       }
                       swaps[j] = maxindex;
                                   
     
                       if(j < (n - 1))
                       {
                            min = A[(j + 1) + (j + 1) * n] - 
                                  A[j + (j + 1) * n] * A[j + (j + 1) * n] / 
                                  A[j + j * n];         
                             
                            for(i = (j + 2); i < n; i++)
                            {
                                  temp = A[i + i * n] - 
                                         A[j + i * n] * A[j + i * n] /
                                         A[j + j * n];
                                  if(temp < min) min = temp;
                            }
                       }

                  if(min < - mu * gamma)
                       {
                              phaseone = 0;
                       }
                       else //j-th iteration of factorization
                       {
                            jthFac(A, L, j, n);
                            j = j + 1;
                       }
               }
       }
       
       //phasetwo, A not positive definite
       double delta;
       double sumI, sumN;
       double delta_old = 0; 
       int kv = j; 

       if(phaseone == 0 && j == (n - 1)) 
       {
             if(- tau * A[n * n - 1] / (1 - tau) > tau_bar * gamma)
             {
                  

                  delta = - A[n * n - 1] - tau * A[n * n - 1] / (1 - tau);
             }
             else
             {
                  delta = - A[n * n - 1] + tau_bar * gamma;
             }
             A[n * n - 1] += delta;
             L[n * n - 1] = sqrt(A[n * n - 1]);
             swaps[n - 1] = n - 1;
       }  
          
       if(phaseone == 0 && j < (n - 1))
       {
             double g[(n - j)];        
             // calculate lower Gerschgorin bounds of A_{k+1}
             for(i = j; i < n; i++)
             {
                   sumI = 0;
                   sumN = 0;
                   
                   for(k = j; k < i; k++)
                   {
                         sumI += fabs(A[i + k * n]);
                   }
                   for(k = i + 1; k < n; k++)          
                   {
                         sumN += fabs(A[k + i * n]);                    
                   }
                   
                   g[i - j] = A[i + i * n] - sumI - sumN;   
             }
             
             
             // Modified Cholesky Decomposition
             for( ; j < (n - 2); j++)
             {
                   // Pivot on maximum lower Gerschgorin bound estimate
                   maxindex = j;        
                   max = g[j - kv];          
                   for(i = j + 1; i < n; i++)
                   {
                         if(g[i - kv] > max)
                         {
                              max = g[i - kv];   
                              maxindex = i;
                         }
                   }
                   
                   if(maxindex != j) 
                   {
                          rowColSwap(A, j, maxindex, n);
                          rowColSwap(L, j, maxindex, n);
                          
                          //swap Gerschgorin bounds;
                          sumI = g[j - kv];
                          g[j - kv] = g[maxindex - kv];
                          g[maxindex - kv] = sumI;
                   }
                   swaps[j] = maxindex;
                   
                   // Calculate E_jj and add to diagonal
                   sumN = 0;
                   for(i = (j + 1); i < n; i++)
                   {
                         sumN += fabs(A[i + j * n]);
                   }
                   
                   delta = 0;
                   if(-A[j + j * n] + sumN > delta) 
                   {
                           delta = -A[j + j * n] + sumN;
                   }
                   if(-A[j + j * n] + tau_bar * gamma > delta)
                   {
                           delta = -A[j + j * n] + tau_bar * gamma;
                   }
                   if(delta_old > delta)
                   {
                           delta = delta_old;
                   }
                   
                   if(delta > 0)
                   {
                            A[j + j * n] += delta;
                            delta_old = delta;
                   }
                   
                   // Update Gerschgorin bound estimates
                   if(A[j + j * n] != sumN)
                   {
                          temp = 1 - sumN / A[j + j * n];
                          for(i = j + 1; i < n; i++)
                          {
                                g[i - k] += fabs(A[i + j * n]) * temp;      
                          }
                   }
                   
                   // Perform jth iteration of factorization
                   jthFac(A, L, j, n);
             }
             
             // Final 2 x 2 submatrix
       
             double lambda_lo = 0;
             double lambda_hi = 0;
       
             // EIGENWERTPROBLEM
             double a = (A[(n - 1) * n - 2] + A[n * n - 1]) / 2;
             double b = sqrt((A[(n - 1) * n - 2] - A[n * n - 1]) * 
                             (A[(n - 1) * n - 2] - A[n * n - 1]) / 4 +
                              A[(n - 1) * n - 1] * A[(n - 1) * n - 1]);
       
             lambda_hi = a + b;
             lambda_lo = a - b;       
       
             delta = 0;
             if(-lambda_lo + tau * (lambda_hi - lambda_lo) / (1 - tau) > delta)
             {
                           delta = -lambda_lo + tau * 
                                   (lambda_hi - lambda_lo) / (1 - tau);         
             }
             if(-lambda_lo + tau_bar * gamma > delta)
             {
                           delta = -lambda_lo + tau_bar * gamma;      
             }
             if(delta_old > delta)
             {
                          delta = delta_old;        
             }
       
             if(delta > 0)
             {
                      A[(n - 1) * n - 2] += delta;
                      A[n * n - 1] += delta;
                      delta_old = delta; 

             }
       
             L[(n - 1) * n - 2] = sqrt(A[(n - 1) * n - 2]); 
             L[n * n - 2] = A[n * n - 2] / L[(n - 1) * n - 2];
             L[n * n - 1] = sqrt(A[n * n - 1] - L[n * n - 2] * L[n * n - 2]);
             // Overwrites A in all 3 cases  
             
             swaps[n - 2] = n - 2;
             swaps[n - 1] = n - 1;
       }
}


SEXP cholesky(SEXP X, SEXP N, SEXP TAU, SEXP TAU_BAR, SEXP MU)
{
     SEXP A = duplicate(X);
     PROTECT(A);
     
     int n = *REAL(N);
     int i;
     double tau = *REAL(TAU);
     double tau_bar = *REAL(TAU_BAR);
     double mu = *REAL(MU);
     double swaps[n];
     
     SEXP L; 
     PROTECT(L = allocVector(REALSXP, n * (n + 1)));
     
     double *l = REAL(L);
     RMCDA(REAL(A), l, n, tau, tau_bar, mu, swaps);
     
     for(i = 0; i < n; i++)
     {
           l[i + n * n] = swaps[i];
     }
     
     UNPROTECT(2);
     return L;
}


SEXP teststat(SEXP Y, SEXP SIGMA)
{
     PROTECT(Y);
     double sigma = *REAL(SIGMA);
     
     SEXP MAX;
     PROTECT(MAX = allocVector(REALSXP, 1));
     double *max = REAL(MAX);
     max[0] = 0;
     
     int n = length(Y);
     
     double sqn = sqrt(n);
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


SEXP teststat_ma(SEXP Y, SEXP SIGMA, SEXP SWAPS, SEXP N, SEXP M)
{
     PROTECT(Y);
     PROTECT(SIGMA);
     PROTECT(SWAPS);
     double *sigma = REAL(SIGMA);
     double *swaps = REAL(SWAPS);
     
     int n = *REAL(N);
     int m = *REAL(M);
     
     /*
      * cumsum: matrix; each column is the cumulative sum of the corresponding 
      * column in Y
      */
     double *cumsum = REAL(cumsum_ma(Y, N, M));
     
     double sumN[m];
     double temp[m];
     double maxCand, temp2;
     
     SEXP MAX; 
     PROTECT(MAX = allocVector(REALSXP, 1));
     double *max = REAL(MAX);
     max[0] = 0;
     
     int i, j, k, index;
     
     for(j = 0; j < m; j++)
     {
           sumN[j] = cumsum[(j + 1) * n - 1] / n;
     }
     
     for(i = 0; i < n; i++)
     {
           for(j = 0; j < m; j++)
           {
                 temp[j] = cumsum[j * n + i] - (i + 1) * sumN[j];
           }
           
           // swap temp according to swaps
           for(j = (m - 1); j >= 0; j--)
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
        if(x[i] < 1) 
    {    
            t = exp(- M_PI * M_PI / (8 * x[i] * x[i]));
            z = 1 + 8 * pow(t, k_max);
            s = M_SQRT_PI * M_SQRT2 / x[i];
            
            for(k = (k_max - 1); k >= 1; k--)
            {
                           z = (1 + pow(t, 8 * k) * z);
            }
            
            x[i] = s * t * z;
        }
        else 
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


SEXP pKSdist(SEXP statistic, SEXP stol)
{
    int n = LENGTH(statistic);
    double tol = asReal(stol);
    SEXP ans = duplicate(statistic);
    KSdist(n, REAL(ans), tol);
    return ans;
}




/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */

#define ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }

double quick_select(double arr[], int k, int n) 
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0; 
    high = n-1; 
    median = k;
    
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}


// compare function
int cmpfun(const void *a, const void *b) 
{
   return ( *(double*)a - *(double*)b );
}


/* kthPair: Finds the k-th biggest element of the set X + Y 
* x and y must be in descending order
* 1 <= k <= n * m
*/
double kthPair(double *x, double *y, int n, int m, int k)
{               
       int i, j, l, sum1, sum2;
       int ji[n];
       
       int L = 0;
       int R = n * m;
       
       int Lb[n]; 
       int Rb[n], P[n], Q[n], w[n];
       double A[n];
       
       double am;
       
       for(i = 0; i < n; i++)
       {
             Lb[i] = 0;
             Rb[i] = m - 1;
       }
       
       while(R - L > n)
       {
               sum1 = 0;
               
               for(i = 0; i < n; i++)
               {
                     if(Lb[i] <= Rb[i])
                     {
                              ji[i] = (Lb[i] + Rb[i]) / 2;
                              A[i] = x[i] + y[ji[i]];
                              w[i] = Rb[i] - Lb[i] + 1;
                              sum1 += w[i];
                     } else
                     {
                           w[i] = 0;
                           A[i] = 0; 
                     }
               }
                                
               qsort(A, n, sizeof(double), cmpfun);
                     
               l = 0; i = 0;
               int median = (sum1 - 1) / 2;
                     
               while(l <= median)
               {
                       l += w[i];
                       i++;
               }
                    
               // weighted median of A with respect to w 
               am = A[i - 1];
               
               
               for(i = 0; i < n; i++)
               {
                     if(x[i] + y[0] <= am) 
                     {
                             P[i] = -1;
                     } else
                     {
                           j = 1;
                           while(j < m && x[i] + y[j] > am) j++;
                           P[i] = j - 1;
                     }
                     
                     if(x[i] + y[m - 1] >= am) 
                     {
                             Q[i] = m;
                     } else
                     {
                           j = m - 2;
                           while(j >= 0 && x[i] + y[j] < am) j--;
                           Q[i] = j + 1;
                     }
               }
               
               sum1 = 0;
               sum2 = 0;
               
               for(i = 0; i < n; i++)
               {
                     sum1 += P[i];
                     sum2 += Q[i];
               }
               sum1 += n;
               
               if(k <= sum1)
               {
                    for(i = 0; i < n; i++)
                    {
                          Rb[i] = P[i];
                    }
               } else if(k > sum2)
               {
                      for(i = 0; i < n; i++)
                      {
                            Lb[i] = Q[i];
                      }
               } else 
               {
                      return(am);
               }
               
               L = 0;
               R = 0;
               
               for(i = 0; i < n; i++)
               {
                     L += Lb[i];
                     R += Rb[i];
               }
               R += n;
       }
       
       /* create new array wA containing all elements that have not been 
        * excluded from possibly being the k-th largest
        */
       l = 0;
       for(i = 0; i < n; i++)
       {
             if(Lb[i] <= Rb[i])
             {
                      l += Rb[i] - Lb[i] + 1;
             }
             
       }
       
       double wA[l];
       j = 0;
       int ii;
       for(i = 0; i < n; i++)
       {
             for(ii = Lb[i]; ii <= Rb[i]; ii++)
             {
                      wA[j] = x[i] + y[ii];
                      j++;
             }
       }
       
       // return the (k - L)-th largest, i.e. the (l - k - L)-th smallest
       return quick_select(wA, l - k + L, l);
}


SEXP KthPair(SEXP X, SEXP Y, SEXP N, SEXP M, SEXP K)
{
     PROTECT(X);
     PROTECT(Y);
     
     int n = *REAL(N);
     int m = *REAL(M);
     int k = *REAL(K);
     
     double *x = REAL(X);
     double *y = REAL(Y);
     
     SEXP RES;
     PROTECT(RES = allocVector(REALSXP, 1));
     double *res = REAL(RES);
     
     res[0] = kthPair(x, y, n, m, k);
     
     UNPROTECT(3);
     
     return RES;
}


