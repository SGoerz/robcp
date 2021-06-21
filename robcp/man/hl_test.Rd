\name{hl_test}
\alias{hl_test}
\title{Hodges-Lehmann Test for Change Points}

\usage{
hl_test(x, method = "subsampling", control = list(), tol = 1e-8)
}

\description{
Performs the Hodges-Lehmann change point test.
}

\arguments{
\item{x}{time series (numeric or ts vector).}
\item{method}{method for estimating the long run variance.}
\item{control}{a list of control parameters (cf. \code{\link{lrv}}).}
\item{tol}{tolerance of the distribution function (numeric), which is used do compute p-values.}
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
...missing...
default values for b_n and b2??
}

\author{
Sheila \enc{Görz}{Goerz}
}

\references{
Dehling, H., Fried, R., and Wendler, M. "A robust method for shift detection in time series." Biometrika 107.3 (2020): 647-660.
}

\seealso{
\code{\link{HodgesLehmann}}, \code{\link{medianDiff}}, \code{\link{lrv}}
}

\examples{
#time series with a structural break at t = 20
Z <- c(rnorm(20, 0), rnorm(20, 2))

hl_test(Z, control = list(overlapping = TRUE, b_n = 5, b2 = 0.05))
}