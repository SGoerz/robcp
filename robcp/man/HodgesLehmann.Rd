\name{HodgesLehmann}
\alias{HodgesLehmann}
\alias{u_hat}
\title{Hodges Lehmann Test Statistic}
\description{
Computes the test statistic for the Hodges-Lehmann change point test.
}
\usage{
HodgesLehmann(x, b_u = "nrd0", method = "kernel", control = list())
u_hat(x, b_u = "nrd0")
}
\arguments{
\item{x}{time series (numeric or \code{ts} vector).}
\item{b_u}{bandwidth for \code{\link{u_hat}}. Either a numeric value or the name
           of a bandwidth selection function (c.f. \code{\link{bw.nrd0}}).}
\item{method}{method of long run variance estimation. Options are \code{"kernel"}, \code{"subsampling"} and \code{"bootstrap"}.}
\item{control}{a list of control parameters for the estimation of the long run variance (cf. \code{\link{lrv}}).}
}

\details{
Let \eqn{n} be the length of the time series. The Hodges-Lehmann test statistic is then computed as
\deqn{\frac{\sqrt{n}}{\hat{\sigma}_n} \max_{1 \leq k < n} \hat{u}_{k,n}(0) \frac{k}{n} \left(1 - \frac{k}{n}\right) | med\{(x_j - x_i); 1 \leq i \leq k; k+1 \leq j \leq n\} | ,}
where \eqn{\hat{\sigma}} is the estimated long run variance, computed by the square root of \code{\link{lrv}}. By default the long run variance is estimated kernel-based with the following bandwidth: first the time series \eqn{x} is corrected for the estimated change point and Spearman's autocorrelation to lag 1 (\eqn{\rho}) is computed. Then the default block length follows as 
\deqn{l = \max\left\{\left\lceil n^{0.45} \left( \frac{2\rho}{1 - \rho^2}\right)^{0.9} \right\rceil, 1\right\}.}

\eqn{\hat{u}_{k,n}(0)} is estimated by \code{u_hat} on data \eqn{\tilde{x}}, where \eqn{med\{(x_j - x_i); 1 \leq i \leq k; k+1 \leq j \leq n\}} was subtracted from \eqn{x_{k+1}, ..., x_n}. Then \code{\link{density}} with the arguments \code{na.rm = TRUE}, \code{from = 0}, \code{to = 0}, \code{n = 1} and \code{bw =  b_u} is applied to \eqn{(\tilde{x}_i - \tilde{x_j})_{1 \leq i < j \leq n}}.
}
%\eqn{\hat{u}_{k,n}(0)} is estimated by \code{u_hat} on data \eqn{\tilde{x}}, where \eqn{med\{(x_j - x_i); 1 \leq i \leq k; k+1 \leq j \leq n\}} was subtracted from \eqn{x_{k+1}, ..., x_n}. Then \code{u_hat} returns
%\deqn{\frac{2}{n(n-1)b_u} \sum_{1 \leq i < j \leq n} QS((\tilde{x}_i - \tilde{x}_j) / b_u),}
%where \eqn{QS} is the Quadratic Spectral kernel function (cf. Section 'Note' in \code{\link{lrv}}). By default, \eqn{b_u} is computed from the set of all possible differences between two values in \code{x} using \code{\link{bw.SJ}} or \code{\link{bw.nrd0}}} 

\value{
\code{HodgesLehmann} returns a test statistic (numeric value) with the following attributes:
\item{cp-location}{indicating at which index a change point is most likely.}
\item{teststat}{test process (before taking the maximum).}
\item{lrv-estimation}{long run variance estimation method.}
\item{sigma}{estimated long run variance.}
\item{param}{parameter used for the lrv estimation.}
\item{kFun}{kernel function used for the lrv estimation.}
Is an S3 object of the class "cpStat".

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
\donttest{HodgesLehmann(x, b_u = 0.01)}
}
