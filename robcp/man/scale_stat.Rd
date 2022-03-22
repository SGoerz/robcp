\name{scale_stat}
\alias{scale_stat}

\title{Test statistic to detect Scale Changes}

\description{
Computes the test statistic for CUSUM-based tests on scale changes.
}

\usage{
scale_stat(x, version =  c("empVar", "MD", "GMD"), method = "kernel",
           control = list(), constant = 1.4826)%, alpha = 0.5)
}

\arguments{
\item{x}{time series (numeric or ts vector).}

\item{version}{variance estimation method.}

\item{method}{either \code{"kernel"} for performing a kernel-based long run variance estimation, or \code{"bootstrap"} for performing a dependent wild bootstrap. See 'Details' below.}

\item{control}{a list of control parameters.}

\item{constant}{scale factor for the MAD.}

%\item{alpha}{quantile of the distribution function of all absolute pairwise differences used in \code{version = "Qalpha"}.}
}

\value{
Test statistic (numeric value) with the following attributes:
\item{cp-location}{indicating at which index a change point is most likely.}
\item{teststat}{test process (before taking the maximum).}
\item{lrv-estimation}{long run variance estimation method.}
If \code{method = "kernel"} the following attributes are also included:
\item{sigma}{estimated long run variance.}
\item{param}{parameter used for the lrv estimation.}
\item{kFun}{kernel function used for the lrv estimation.}

Is an S3 object of the class "cpStat".
}

\details{
Let \eqn{n} be the length of the time series. The CUSUM test statistic for testing on a change in the scale is then defined as
\deqn{\hat{T}_{s} = \max_{1 < k \leq n} \frac{k}{\sqrt{n}} |\hat{s}_{1:k} - \hat{s}_{1:n}|,}
where \eqn{\hat{s}_{1:k}} is a scale estimator computed using only the first \eqn{k} elements of the time series \eqn{x}.

If \code{method = "kernel"}, the test statistic \eqn{\hat{T}_s} is divided by the estimated long run variance \eqn{\hat{D}_s} so that it asymptotically follows a Kolmogorov distribution. \eqn{\hat{D}_s} is computed by the function \code{\link{lrv}} using kernel-based estimation. 

For the scale estimator \eqn{\hat{s}_{1:k}}, there are five different options which can be set via the \code{version} parameter:

\bold{Empirical variance} (\code{empVar})
\deqn{\hat{\sigma}^2_{1:k} = \frac{1}{k-1} \sum_{i = 1}^k (x_i - \bar{x}_{1:k})^2; \; \bar{x}_{1:k} = k^{-1} \sum_{i = 1}^k x_i.}

\bold{Mean deviation} (\code{MD})
\deqn{\hat{d}_{1:k}= \frac{1}{k-1} \sum_{i = 1}^k |x_i - med_{1:k}|,}
where \eqn{med_{1:k}} is the median of \eqn{x_1, ..., x_k}. 

\bold{Gini's mean difference} (\code{GMD})
\deqn{\hat{g}_{1:k} = \frac{2}{k(k-1)} \sum_{1 \leq i < j \leq k} (x_i - x_j).}

%\bold{Median absolute deviation} (\code{MAD})
%\deqn{\hat{m}_{1:k} = med(x_1 - med_{1:k}, ..., x_k - med_{1:k}).}

%\bold{\eqn{Q^{\alpha}}} (\code{Qalpha})
%\deqn{\hat{Q}^{\alpha}_{1:k} = U_{1:k}^{-1}(\alpha) = \inf\{x | \alpha \leq U_{1:k}(x)\},}
%where \eqn{U_{1:k}} is the empirical distribtion function of \eqn{|x_i - x_j|, \, 1 \leq i < j \leq k} (cf. \code{\link{Qalpha}}). 


For the kernel-based long run variance estimation, the default bandwidth \eqn{b_n} is determined as follows:

If \eqn{\hat{\rho}_j} is the estimated autocorrelation to lag \eqn{j}, a maximal lag \eqn{l} is selected to be the smallest integer \eqn{k} so that 
\deqn{\max \{|\hat{\rho}_k|, ..., |\hat{\rho}_{k + \kappa_n}|\} \leq 2 \sqrt(\log_{10}(n) / n), }
\eqn{\kappa_n = \max \{5, \sqrt{\log_{10}(n)}\}}. This \eqn{l} is determined for both the original data \eqn{x} and the squared data \eqn{x^2} and the maximum \eqn{l_{max}} is taken. Then the bandwidth \eqn{b_n} is the minimum of \eqn{l_{max}} and \eqn{n^{1/3}}.

}

\author{
Sheila \enc{Görz}{Goerz}
}

\references{
Gerstenberger, C., Vogel, D., and Wendler, M. (2020). Tests for scale changes based on pairwise differences. Journal of the American Statistical Association, 115(531), 1336-1348.
}

\seealso{
\code{\link{lrv}}, \code{\link{Qalpha}}
}

\examples{
x <- arima.sim(list(ar = 0.5), 100)

# under H_0:
scale_stat(x, "GMD")
%scale_stat(x, "Qalpha", method = "bootstrap")

# under the alternative:
x[51:100] <- x[51:100] * 3
scale_stat(x)
%scale_stat(x, "Qalpha", method = "bootstrap")
}
