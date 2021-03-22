\name{lrv}
\alias{lrv}

\title{
Long Run Variance
}
\description{
Estimates the long run variance respectively covariance matrix of the supplied time series via kernel estimation with a flat top kernel.
}

\usage{
lrv(x, b_n, gamma0 = TRUE)
}

\arguments{
  \item{x}{vector or matrix with each column representing a time series (numeric).}
  \item{b_n}{Bandwidth for kernel-based computation. How many estimated autocovariances should be used to compute the long run variance estimation? Must be greater than 0. Default is \eqn{n^{1/3}} with \eqn{n} being the number of observations.} 
  \item{gamma0}{boolean: should only the variance of the time series, i.e. \eqn{\gamma}(0), be returned if the long run variance was estimated to be < 0?}
}

\details{
The long run variance equals \eqn{n} times the asymptotic variance of the arithmetic mean of a short range dependent time series, where \eqn{n} is the length of the time series. It is used to standardize CUSUM Tests.

The long run variance is estimated by a kernel estimator using the bandwidth \eqn{b_n = n^{1/3}} and the flat top kernel \deqn{k(x) = 1 * 1_{|x| <= 0.5} + (2 - 2 * |x|) * 1_{0.5 < |x| < 1}.} In the one dimensional case this results in:

\deqn{\hat{\sigma}^2 = (1 / n) \sum_{i = 1}^n (x[i] - mean(x))^2 + (2 / n) \sum_{h = 1}^{n - 1} \sum_{i = 1}^{n - h} (x[i] - mean(x)) * (x[i + h] - mean(x)) * k(h / b_n).}

If \code{x} is a multivariate timeseries the \eqn{k,l}-element of \eqn{\Sigma} is estimated by

\deqn{\hat{\Sigma}^{(k,l)} = (1 / n) \sum_{i,j = 1}^{n}(x[i]^{(k)} - mean(x)^{(k)}) * (x[j]^{(l)} - mean(x)^{(l)}) * k((i-j) / b_n).}
}

\value{
long run variance \eqn{\sigma^2} resp. \eqn{\Sigma} (numeric)
}

\author{
Sheila \enc{Görz}{Goerz}
}

\seealso{
\code{\link{psi}},
\code{\link{psi_cumsum}}, \code{\link{CUSUM}}, \code{\link{pKSdist}}, 
\code{\link{huber_cusum}}
}

\examples{
Z <- c(rnorm(20), rnorm(20, 2))
sigma2(Z)
}
