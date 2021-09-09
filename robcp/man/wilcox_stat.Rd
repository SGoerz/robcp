\name{wilcox_stat}
\alias{wilcox_stat}
\title{Wilcoxon-Mann-Whitney Test Statistic for Change Points}
\description{
Computes the test statistic for the Wilcoxon-Mann-Whitney change point test
}

\usage{
wilcox_stat(x, h = 1L, method = "subsampling", control = list())
}

\arguments{
\item{x}{Time series (numeric or ts vector)}
\item{h}{Kernel function of the U statistic (1L or 2L, or a function with two parameters).}
\item{method}{Method for estimating the long run variance}
\item{control}{A list of control parameters for the estimation of the long run variance (cf. \code{\link{lrv}})}
}

\details{
Let n be the length of \code{x}, i.e. the number of observations.

\code{h = 1L}:
\deqn{T_n = \frac{1}{\hat{\sigma}} \max_{1 \leq k \leq n} \left| \frac{1}{n^{3/2}} \sum_{i = 1}^k \sum_{j = k+1}^n (1_{\{x_i < x_j\}} - 0.5) \right|}


\code{h = 2L}:
\deqn{T_n = \frac{1}{\hat{\sigma}} \max_{1 \leq k \leq n} \left| \frac{1}{n^{3/2}} \sum_{i = 1}^k \sum_{j = k+1}^n (x_i - x_j) \right|}

\eqn{\hat{\sigma}} is estimated by the square root of \code{\link{lrv}}. The denominator corresponds to that in the ordinary CUSUM change point test. 

By default, non-overlapping subsampling is used. If \code{h = 1L}, the default for \code{distr} is \code{TRUE} and otherwise it is \code{FALSE}. If no block length is supplied, first the time series \eqn{x} is corrected for the estimated change point and Spearman's autocorrelation to lag 1 (\eqn{\rho}) is computed. Then the default block length follows as 
\deqn{l = \max\left\{\left\lceil n^{1/3} \left( \frac{2\rho}{1 - \rho^2}\right)^{2/3} \right\rceil, 1\right\}.}
}

\value{
Test statistic (numeric value) with the attribute cp-location indicating at which index a change point is most likely. Is an S3 object of the class \code{cpStat}
}

\references{
Dehling, H., et al. "Change-point detection under dependence based on two-sample U-statistics." Asymptotic laws and methods in stochastics. Springer, New York, NY, 2015. 195-220.
}

\seealso{
\code{\link{lrv}}
}

\author{
Sheila \enc{Görz}{Goerz}
}

\examples{
# time series with a location change at t = 20
x <- c(rnorm(20, 0), rnorm(20, 2))

# Wilcoxon-Mann-Whitney change point test statistic
wilcox_stat(x, h = 1L, control = list(b_n = length(x)^(1/3)))
}