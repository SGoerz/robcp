\name{hl_test}
\alias{hl_test}
\title{Hodges-Lehmann Test for Change Points}

\usage{
hl_test(x, b_u = "SJ", method = "subsampling", control = list(), tol = 1e-8)
}

\description{
Performs the two-sample Hodges-Lehmann change point test.
}

\arguments{
\item{x}{time series (numeric or ts vector).}
\item{b_u}{bandwidth for \code{\link{u_hat}}.}
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
The function performs the two-sample Hodges-Lehmann change point test. It tests the hypothesis pair
\deqn{H_0: \mu_1 = ... = \mu_n}
\deqn{vs.}
\deqn{H_1: \exists k \in \{1, ..., n-1\}: \mu_k \neq \mu_{k+1}}
where \eqn{\mu_t = E(X_t)} and \eqn{n} is the length of the time series. \eqn{k} is called a 'change point'. 

The test statistic is computed using \code{\link{HodgesLehmann}} and asymptotically follows a Kolmogorov distribution. To derive the p-value, the function \code{\link{pKSdist}} is used. 
}

\author{
Sheila \enc{Görz}{Goerz}
}

\references{
Dehling, H., Fried, R., and Wendler, M. "A robust method for shift detection in time series." Biometrika 107.3 (2020): 647-660.
}

\seealso{
\code{\link{HodgesLehmann}}, \code{\link{medianDiff}}, \code{\link{lrv}}, \code{\link{pKSdist}}
}

\examples{
#time series with a structural break at t = 20
Z <- c(rnorm(20, 0), rnorm(20, 2))

hl_test(Z, control = list(overlapping = TRUE, b_n = 5, b2 = 0.05))
}