\encoding{utf8}
\name{CUSUM}
\alias{CUSUM}

\title{
CUSUM Test Statistic
}

\description{
Computes the test statistic for the CUSUM change point test.
}

\usage{
CUSUM(x, method = "kernel", control = list(), inverse = "Cholesky", ...)
}

\arguments{
	\item{x}{vector or matrix with each column representing a time series (numeric).}
	\item{method}{method of long run variance estimation. Options are \code{"kernel"}, \code{"subsampling"} and \code{"bootstrap"}.}
	\item{control}{a list of control parameters for the estimation of the long run variance (cf. \code{\link{lrv}}).}
	\item{inverse}{character string specifying the method of inversion. Options are "Cholesky" for inverting over \code{\link{modifChol}} and "generalized" for using \code{\link[MASS]{ginv}} from the \code{MASS} package.}
	\item{...}{further arguments passed to the inverse-computing functions.}
}

\value{
Test statistic (numeric value) with the attribute \code{cp-location} indicating at which index a change point is most likely.
}

\details{
Let n be the length of the time series \eqn{x = (x_1, ..., x_n)}. 

In case of a univariate time series the test statistic can be written as \deqn{\max_{k = 1, ..., n}\frac{1}{\sqrt{n} \sigma}\left|\sum_{i = 1}^{k} x_i - (k / n) \sum_{i = 1}^n x_i\right|,} where \eqn{\sigma} is the square root of \code{\link{lrv}}.
Default method is \code{"kernel"} and the default kernel function is \code{"TH"}. If no bandwidth value is supplied, first the time series \eqn{x} is corrected for the estimated change point and Spearman's autocorrelation to lag 1 (\eqn{\rho}) is computed. Then the default bandwidth follows as 
\deqn{b_n = \max\left\{\left\lceil n^{1/3} \left( \frac{2\rho}{1 - \rho^2} \right)^{2/3} \right\rceil, 1 \right\}.}


In case of a multivariate time series the test statistic follows as
\deqn{\max_{k = 1, ..., n}\frac{1}{n}\left(\sum_{i = 1}^{k} X_i - \frac{k}{n} \sum_{i = 1}^{n} X_i\right)^T \Sigma^{-1} \left(\sum_{i = 1}^{k} X_i - \frac{k}{n} \sum_{i = 1}^{n} X_i\right),} where \eqn{X_i} denotes the i-th row of x and \eqn{\Sigma^{-1}} is the inverse of \code{\link{lrv}}.
}

\seealso{
\code{\link{psi_cumsum}}, 
\code{\link{psi}}
}


\author{
Sheila \enc{Görz}{Goerz}
}


\examples{
# time series with a location change at t = 20
ts <- c(rnorm(20, 0), rnorm(20, 2))

# Huberized CUSUM change point test statistic
CUSUM(psi(ts))
}
