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
\item{h}{version of the test (integer, 1 or 2)}
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

}

\author{
Sheila \enc{Görz}{Goerz}
}

\references{
Dehling, H., et al. "Change-point detection under dependence based on two-sample U-statistics." Asymptotic laws and methods in stochastics. Springer, New York, NY, 2015. 195-220.
}

\seealso{
\code{\link{wilcox_stat}}, \code{\link{lrv}}
}

\examples{
#time series with a structural break at t = 20
Z <- c(rnorm(20, 0), rnorm(20, 2))

wmw_test(Z, h = 1L, control = list(overlapping = TRUE, b_n = 5))
}
