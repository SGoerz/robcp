#include <Rcpp.h>
#include <vector>
#include <R_ext/Arith.h>
using namespace Rcpp;

bool sortcol(const std::vector<double>& v1, const std::vector<double>& v2)
{
  return v1[0] > v2[0];
}

//' Computes the weighted median of a numeric vector.
//' @param x Numeric vector.
//' @param w Integer vector of weights
//'
//' @return Weighted median of x with respect to w.
//'
//' @examples 
//' x <- c(1, 4, 9)
//' w <- c(5, 1, 1)
//' weightedMedian(x, w)
// [[Rcpp::export]]
double weightedMedian(NumericVector x, IntegerVector w)
{
  int n = w.size();
  
  if(n != x.size()) Rcpp::stop("x and w need to have the same length!");
  int i, l;
  
  bool allPos = true;
  for(i = 0; i < n; i++)
  {
    if(w[i] < 0) allPos = false;
  }
  if(!allPos) Rcpp::stop("Negative weights supplied!");
  
  int med = 0;
  
  std::vector<std::vector<double>> vect(n, std::vector<double> (2, 0));
  
  for(i = 0; i < n; i++)
  {
    vect[i][0] = x[i];
    vect[i][1] = w[i];
    med += w[i];
  }
  
  std::sort(vect.begin(), vect.end(), sortcol);
  
  l = 0; i = 0;
  med = (med + 1) / 2;
  
  while(l < med)
  {
    l += vect[i][1];
    i++;
  }
  
  // weighted median of A with respect to w
  return vect[i - 1][0];
}

//'kthPair:
//'
//'input: - x and y (numeric vectors; descending order)
//'       - n and m (lengths of x and y; integer)
//'       - k (index of element to choose; integer; 1 <= k <= n * m)
//'
// [[Rcpp::export]]
double kthPair(NumericVector x1, NumericVector y1, int k, int k2 = NA_INTEGER)
{
  NumericVector x = clone(x1);
  NumericVector y = clone(y1);
  
  int n = x.size();
  int m = y.size();
  
  if(k <= 0 || k > n * m) Rcpp::stop("k out of bounds");
  if(IntegerVector::is_na(k2))
  {
    k2 = k;
  } else if(k2 <= 0 || k2 > n * m)
  {
    stop("k2 out of bounds");
  } else if(std::abs(k - k2) > 1)
  {
    Rcpp::stop("k and k2 must be consecutive indices!");
  }

  NumericVector temp(2);
  temp[0] = R_NaReal;
  temp[1] = R_NaReal;

  std::sort(x.begin(), x.end());
  std::sort(y.begin(), y.end());
  
  std::reverse(x.begin(), x.end());
  std::reverse(y.begin(), y.end());
  
  int i, j, l, sum1, sum2;
  IntegerVector ji(n);
  
  int L = 0;
  int R = n * m;
  
  IntegerVector Lb(n), Rb(n), P(n), Q(n), w(n);
  NumericVector A(n);
  
  double am;
  
  for(i = 0; i < n; i++)
  {
    Lb[i] = 0;
    Rb[i] = m - 1;
  }
  
  while(R - L > n)
  {
    Rcpp::checkUserInterrupt();
    
    sum1 = 0;
    
    // find 'middle' elements per row of possible candidates
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
    
    am = weightedMedian(A, w);

    // part values in sets higher than median; containing median; 
    // and lower than median
    // P: border after values; Q: border before value
    
    j = std::min(n-1, Lb[n-1]);
    
    for(i = n-1; i >= 0; i--)
    {
      while(j < m && x[i] + y[j] > am) j++;
      P[i] = j - 1;
    }  

    j = std::max(0, Rb[0]);
    for(i = 0; i < n; i++)
    {
      while(j >= 0 && x[i]  + y[j] < am) j--;
      Q[i] = j + 1;
    }
    
    sum1 = 0;
    sum2 = 0;
    
    // P contains sum1 values; Q excludes sum2 values
    for(i = 0; i < n; i++)
    {
      sum1 += P[i];
      sum2 += Q[i];
    }
    sum1 += n;
    
    // check in which set to search
    if(k <= sum1 || k2 <= sum1)
    {
      for(i = 0; i < n; i++)
      {
        Rb[i] = P[i];
      }
    } else if(k > sum2 || k2 > sum2)
    {
      for(i = 0; i < n; i++)
      {
        Lb[i] = Q[i];
      }
    } 
    
    if(k > sum1 && k <= sum2)
    {
      temp[0] = am;
    }
    if(k2 > sum1 && k2 <= sum2)
    {
      temp[1] = am;
    }
    if(!NumericVector::is_na(temp[0]) && !NumericVector::is_na(temp[1]))
    {
      return (temp[0] + temp[1]) / 2;
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
  
  NumericVector wA(l);
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
  
  // return the (k - L)-th largest, i.e. the (l - k + L)-th smallest 
  std::sort(wA.begin(), wA.end());
  
  if(NumericVector::is_na(temp[0]))
  {
    temp[0] = wA[l - k + L];
  }
  if(NumericVector::is_na(temp[1]))
  {
    temp[1] = wA[l - k2 + L];
  }
  
  return (temp[0] + temp[1]) / 2;
  //return wA[l - k + L];
}


