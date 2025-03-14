\name{cor_stat}
\alias{cor_stat}

\title{Test statistic to detect Correlation Changes}

\description{
Computes the test statistic for a CUSUM-based tests on changes in Spearman's rho or Kendall's tau.
}

\usage{
cor_stat(x, version = c("tau", "rho"), method = "kernel", control = list())
}

\arguments{
\item{x}{time series (numeric or ts vector).}
\item{version}{version of the test. Either \code{"rho"} or \code{"tau"}.}
\item{method}{methods of long run variance estimation. Options are \code{"kernel"} and \code{"none"}.}
\item{control}{a list of control parameters.}
}

\value{
Test statistic (numeric value) with the following attributes:
\item{cp-location}{indicating at which index a change point is most likely.}
\item{teststat}{test process (before taking the maximum).}
\item{lrv-estimation}{long run variance estimation method.}
\item{sigma}{estimated long run variance.}
\item{param}{parameter used for the lrv estimation.}
\item{kFun}{kernel function used for the lrv estimation.}

Is an S3 object of the class "cpStat".
}

\details{
Let \eqn{n} be the length of the time series, i.e. the number of rows in \code{x}. In general, the (scaled) CUSUM test statistic is defined as
\deqn{\hat{T}_{\xi; n} = \max_{k = 1, ..., n} \frac{k}{2\sqrt{n}\hat{\sigma}} | \hat{\xi}_k - \hat{\xi}_n |,}
where \eqn{\hat{\xi}} is an estimator for the property on which to test, and \eqn{\hat{\sigma}} is an estimator for the square root of the corresponding long run variance (cf. \code{\link{lrv}}).


If \code{version = "tau"}, the function tests if the correlation between \eqn{x_i} and \eqn{x_i} of the bivariate time series \eqn{(x_i, x_i)_{i = 1, ..., n}} stays constant for all \eqn{i = 1, ..., n} by considering Kendall's tau. Therefore, \eqn{\hat{\xi} = \hat{\tau}} is the the sample version of Kendall's tau:
\deqn{\hat{\tau}_k = \frac{2}{k(k-1)} \sum_{1 \leq i < j \leq k} sign\left((x_j - x_i)(y_j - y_i)\right).}
The default bandwidth for the kernel-based long run variance estimation is \eqn{b_n = \lfloor 2n^{1/3} \rfloor} and the default kernel function is the quatratic kernel.

If \code{version = "rho"}, the function tests if the correlation of a time series of an arbitrary dimension \eqn{d} (>= 2) stays constant by considering a multivariate version of Spearman's rho. Therefore, \eqn{\hat{\xi} = \hat{\rho}} is the sample version of Spearman's rho:
\deqn{\hat{\rho}_k = a(d) \left( \frac{2^d}{k} \sum_{j = 1}^k \prod_{i = 1}^d (1 - U_{i, j; n}) - 1 \right)}
where \eqn{U_{i, j; n} = n^{-1}} (rank of \eqn{x_{i,j}} in \eqn{x_{i,1}, ..., x_{i,n})} and \eqn{a(d) = (d+1) / (2^d - d - 1)}. Here it is essential to use \eqn{\hat{U}_{i, j; n}} instead of \eqn{\hat{U}_{i, j; k}}. The default bandwidth for the kernel-based long run variance estimation is \eqn{\sqrt{n}} and the default kernel function is the Bartlett kernel.

}

\author{
Sheila \enc{Görz}{Goerz}
}

\references{
Wied, D., Dehling, H., Van Kampen, M., and Vogel, D. (2014). A fluctuation test for constant Spearman’s rho with nuisance-free limit distribution. Computational Statistics & Data Analysis, 76, 723-736.
}

\seealso{
\code{\link{lrv}}, \code{\link{cor_cusum}}
}
