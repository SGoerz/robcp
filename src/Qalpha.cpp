#include <Rcpp.h>
#include <R_ext/Arith.h>
using namespace Rcpp;

//' Computes the Q-alpha statistic.
//' @param x Numeric vector.
//' @param alpha quantile. Numeric value in (0, 1]
//'
//' @return numeric vector of Qalpha-s estimated using x[1], ..., x[k], k = 1, ..., n, n being the length of x.
//'
//' @examples 
//' x <- rnorm(10)
//' Qalpha(x, 0.5)
// [[Rcpp::export]]
NumericVector Qalpha(NumericVector x, double alpha = 0.8)
{
  if(alpha <= 0 || alpha > 1) Rcpp::stop("alpha need to be from the interval (0, 1]!");
  
  int n = x.size();
  NumericVector out(n-1);
  NumericVector diffs(n * (n-1) / 2);
  
  int k, i, j;
  
  j = -1;
  
  for(k = 1; k < n; k++)
  {
    for(i = 0; i < k; i++)
    {
      j++;
      diffs[j] = std::abs(x[i] - x[k]);
    }
    
    if(j > 0)
    {
      std::sort(diffs.begin() + j - k + 1, diffs.begin() + j + 1);
      std::inplace_merge(diffs.begin(), diffs.begin() + j - k + 1, 
                         diffs.begin() + j + 1);
    }
    
    out[k-1] = diffs[std::ceil((j+1) * alpha) - 1];
  }
  
  return(out);
} 