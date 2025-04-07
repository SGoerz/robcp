#include <Rcpp.h>
using namespace Rcpp;
// 
// // [[Rcpp::export]]
// NumericVector CUSUM_cpp(NumericVector x)
// {
//   int n = x.size();
//   NumericVector y(n-1);
//   NumericVector csum = cumsum(x);
//   double meanN = csum[n-1] / n;
//   double sqn = sqrt(n);
//   
//   for(int i = 0; i < n-1; i++)
//   {
//     y[i] = fabs(csum[i] - (i + 1) * meanN) / sqn;
//   }
//   
//   return y; 
// }
// 
// // [[Rcpp::export]]
// NumericMatrix cumsum_ma_cpp(NumericMatrix X)
// {
//   int n = X.nrow();
//   int m = X.cols();
//   
//   NumericMatrix Y = clone(X);
//   
//   for(int j = 0; j < m; j++)
//   {
//     for(int i = 1; i < n; i++)
//     {
//       Y(i, j) += Y(i - 1, j);
//     }
//   }
//   
//   return Y;
// }
// 
// // [[Rcpp::export]]
// NumericVector CUSUM_ma_cpp(NumericMatrix X, NumericMatrix Sigma)
// {
//   int n = X.nrow();
//   int m = X.ncol();
//   
//   NumericMatrix csum = cumsum_ma_cpp(X);
//   NumericVector temp(m);
//   NumericVector Y(n);
//   
//   int i, j, k;
//   
//   for(i = 0; i < n; i++)
//   {
//     for(j = 0; j < m; j++)
//     {
//       // temp: CUSUM statistic for each column
//       temp[j] = csum(i, j) - (i + 1) * csum(n-1, j) / n;
//     }
//     
//     Y[i] = 0;
//     
//     for(j = 0; j < m; j++)
//     {
//       for(k = j; k < m; k++)
//       {
//         if(j == k)
//         {
//           Y[i] += temp[j]*temp[j] * Sigma(j, j);
//         }
//         else
//         {
//           Y[i] += 2 * (temp[j]*temp[k] * Sigma(j, k));
//         }
//       }
//     }
//     Y[i] /= n;
//   }
//   
//   return Y;
// }
// 
// // [[Rcpp::export]]
// NumericVector CUSUM_var_cpp(NumericVector X, NumericVector X2)
// {
//   int n = X.size();
//   
//   NumericVector res(n-2);
//   double sqn = sqrt(n);
//   
//   NumericVector csum = cumsum(X);
//   NumericVector csum2 = cumsum(X2);
//   
//   double sumN = pow(csum[n - 1] / n, 2);
//   double sumN2 = csum2[n - 1] / n;
//   
//   int i;
//   
//   for(i = 1; i < n-1; i++)
//   {
//     res[i-1] = fabs(csum2[i] - pow(csum[i], 2) / (i + 1) - (i + 1) * sumN2 + (i + 1) * sumN) / sqn;
//   }
//   
//   return res;
// }
// 
// 
// // [[Rcpp::export]]
// NumericVector MD_cpp(NumericVector X, NumericVector cummed)
// {
//   int n = X.size();
//   NumericVector res(n-1);
//   
//   int i, k;
//   
//   for(k = 1; k < n; k++)
//   {
//     res[k-1] = 0;
//     for(i = 0; i <= k; i++)
//     {
//       res[k-1] += fabs(X[i] - cummed[k]);
//     }
//   }
//   
//   return res;
// }
// 
// // [[Rcpp::export]]
// NumericVector GMD_cpp(NumericVector X)
// {
//   int n = X.size();
//   NumericVector res(n-1);
//   
//   int i, k;
//   
//   res[0] = fabs(X[0] - X[1]);
//   
//   for(k = 2; k < n; k++)
//   {
//     res[k-1] = res[k-2];
//     for(i = 0; i < k; i++)
//     {
//       res[k-1] += fabs(X[i] - X[k]);
//     }
//   }
//   
//   return res;
// }
