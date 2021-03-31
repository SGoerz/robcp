#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>

//** Following: Revised Modified Cholesky Decomposition Algorithm **//

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
// Schnabel, R. B., & Eskow, E. (1999). "A revised modified Cholesky factorization algorithm" 
//                                       SIAM Journal on optimization, 9(4), 1135-1148.
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
    
    // EIGEN VALUE PROBLEM
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
    
    //???
    //swaps[n - 2] = n - 2;
    //swaps[n - 1] = n - 1;
  }
}

// wrapper for the revised modified Cholesky factorization
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
  
  for(i = 0; i < n; i++)
  {
    swaps[i] = i;
  }
  
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