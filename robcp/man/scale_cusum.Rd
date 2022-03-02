\name{scale_cusum}
\alias{scale_cusum}

\title{Tests for Scale Changes Based on Pairwise Differences}

\description{
Performs the CUSUM-based test on changes in the scale.
}

\usage{
scale_cusum(x, version = "empVar", method = "kernel", control = list(), 
            constant = 1.4826, alpha = 0.5, fpc = TRUE, tol = 1e-08)
}

\arguments{
\item{x}{time series (numeric or ts vector).}

\item{version}{variance estimation method (see \code{\link{scale_stat}}).}

\item{method}{method for estimating the long run variance.}

\item{control}{a list of control parameters.}

\item{constant}{scale factor for the MAD.}

\item{alpha}{quantile of the distribution function of all absolute pairwise differences used in \code{version = "Qalpha"}.}

\item{fpc}{finite population correction (boolean).}

\item{tol}{tolerance of the distribution function (numeric), which is used do compute p-values.}
}

\value{
A list of the class "htest" containing the following components:
  \item{statistic}{return value of the function \code{\link{scale_stat}}.}
  \item{p.value}{p-value (numeric).}
  \item{alternative}{alternative hypothesis (character string).}
  \item{method}{name of the performed test (character string).}
  \item{cp.location}{index of the estimated change point location (integer).}
  \item{data.name}{name of the data (character string).}
  
Plus if \code{method = "kernel"}:
\item{lrv}{list containing the compontents \code{method}, \code{param} and \code{value}.}
else if \code{method = "bootstrap"}:
\item{......................}{.................................................................}
}

\details{
The function performs a CUSUM-type test on changes in the scale. Formally, the hypothesis pair is
\deqn{H_0: \sigma^2_2 = ... = \sigma_n^2}
\deqn{vs.}
\deqn{H_1: \exists k \in \{2, ..., n-1\}: \sigma^2_k \neq \sigma^2_{k+1}}
where \eqn{\sigma^2_t = Var(X_t)} and \eqn{n} is the length of the time series. \eqn{k} is called a 'change point'. The hypotheses do not include \eqn{\sigma^2_1} since then the variance of one single observation would need to be estimated. 

The test statistic is computed using \code{\link{scale_stat}} and in the case of \code{method = "kernel"} asymptotically follows a Kolmogorov distribution. To derive the p-value, the function \code{\link{pKSdist}} is used.
}

\references{
Gerstenberger, C., Vogel, D., and Wendler, M. (2020). Tests for scale changes based on pairwise differences. Journal of the American Statistical Association, 115(531), 1336-1348.

\enc{Dürre}{Duerre}, A. (2022+). "Finite sample correction for cusum tests", \emph{unpublished manuscript}
}

\author{
Sheila \enc{Görz}{Goerz}
}

\seealso{
\code{\link{scale_stat}}, \code{\link{lrv}}, \code{\link{pKSdist}}, \code{\link{Qalpha}}
}

\examples{
x <- arima.sim(list(ar = 0.5), 100)

# under H_0:
scale_cusum(x)
scale_cusum(x, "MD")

# under the alternative:
x[51:100] <- x[51:100] * 3
scale_cusum(x)
scale_cusum(x, "MD")
}