\name{HodgesLehmann}
\alias{HodgesLehmann}
\alias{u_hat}
\title{Hodges Lehmann Test Statistic}
\description{
Computes the test statistic for the Hodges-Lehmann change point test.
}
\usage{
HodgesLehmann(x, method = "subsampling", control = list())
u_hat(x, b, kFun = "bartlett")
}
\arguments{
\item{x}{time series (numeric or \code{ts} vector).}
\item{method}{method of long run variance estimation. Options are \code{"kernel"}, \code{"subsampling"} and \code{"bootstrap"}.}
\item{control}{a list of control parameters for the estimation of the long run variance (cf. \code{\link{lrv}}).}
\item{kFun}{character string, name of a kernel function used in \code{u_hat} and for kernel-based estimation in \code{\link{lrv}}.}
\item{b}{bandwidth for \code{u_hat}}
}

\details{
Let \eqn{n} be the length of the time series. The Hodges-Lehmann test statistic is then computed as
\deqn{\frac{\sqrt{n}}{\hat{\sigma}_n} \max_{1 \leq k < n} \hat{u}_{k,n}(x) \frac{k}{n} \left(1 - \frac{k}{n}\right) | med\{(x_j - x_i); 1 \leq i \leq k; k+1 \leq j \leq n\} | ,}
where \eqn{\hat{sigma}} is estimated by the square root of \code{\link{lrv}}. By default the long run variance is estimated via overlapping subsampling.

\eqn{\hat{u}_{k,n}(x)} is estimated by \code{u_hat} on data \eqn{\tilde{x}}, where \eqn{med\{(x_j - x_i); 1 \leq i \leq k; k+1 \leq j \leq n\}} was subtracted from x[(k+1):n]. Then \code{u_hat} returns
\deqn{\frac{2}{n(n-1)b} \sum_{1 \leq i < j \leq n} kFun((\tilde{x}_i - \tilde{x}_j) / b).}
}

\value{
\code{HodgesLehmann} returns a test statistic (numeric value) with the attribute cp-location indicating at which index a change point is most likely. The returned object is an S3 object of the class \code{cpStat}.

\code{u_hat} returns a numeric value.
}

\author{
Sheila \enc{Görz}{Goerz}
}

\references{
Dehling, H., Fried, R., and Wendler, M. "A robust method for shift detection in time series." Biometrika 107.3 (2020): 647-660.
}

\seealso{
\code{\link{medianDiff}}, \code{\link{lrv}}
}

\examples{
# time series with a location change at t = 20
x <- c(rnorm(20, 0), rnorm(20, 2))

# Hodges-Lehmann change point test statistic
HodgesLehmann(x, control = list(b_n = length(x)^(1/3)))
}
