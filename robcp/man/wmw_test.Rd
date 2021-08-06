\name{wmw_test}
\alias{wmw_test}
\title{Wilocxon-Mann-Whitney Test for Change Points}

\usage{
wmw_test(x, h = 1L, method = "subsampling", control = list(), tol = 1e-8)
}

\description{
Performs the Wilcoxon-Mann-Whitney change point test.
}

\arguments{
\item{x}{time series (numeric or ts vector).}
\item{h}{version of the test (integer, 1L or 2L)}
\item{method}{method for estimating the long run variance.}
\item{control}{a list of control parameters (cf. \code{\link{lrv}}).}
\item{tol}{tolerance of the distribution function (numeric), which is used to compute p-values.}
}

\value{
A list of the class "htest" containing the following components:
  \item{statistic}{value of the test statistic (numeric).}
  \item{p.value}{p-value (numeric).}
  \item{alternative}{alternative hypothesis (character string).}
  \item{method}{name of the performed test (character string).}
  \item{cp.location}{index of the estimated change point location (integer).}
  \item{data.name}{name of the data (character string).}
}

\details{
The function performs a Wilcoxon-Mann-Whitney change point test. It tests the hypothesis pair
\deqn{H_0: \mu_1 = ... = \mu_n}
\deqn{vs.}
\deqn{H_1: \exists k \in \{1, ..., n-1\}: \mu_k \neq \mu_{k+1}}
where \eqn{\mu_t = E(X_t)} and \eqn{n} is the length of the time series. \eqn{k} is called a 'change point'. 

The test statistic is computed using \code{\link{wilcox_stat}} and asymptotically follows a Kolmogorov distribution. To derive the p-value, the function \code{\link{pKSdist}} is used. 
}

\author{
Sheila \enc{Görz}{Goerz}
}

\references{
Dehling, H., et al. "Change-point detection under dependence based on two-sample U-statistics." Asymptotic laws and methods in stochastics. Springer, New York, NY, 2015. 195-220.
}

\seealso{
\code{\link{wilcox_stat}}, \code{\link{lrv}}, \code{\link{pKSdist}}
}

\examples{
#time series with a structural break at t = 20
Z <- c(rnorm(20, 0), rnorm(20, 2))

wmw_test(Z, h = 1L, control = list(overlapping = TRUE, b_n = 5))
}
