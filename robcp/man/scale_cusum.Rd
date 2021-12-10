\name{scale_cusum}
\alias{scale_cusum}

\title{Tests for Scale Changes Based on Pairwise Differences}

\description{
Performs the CUSUM-based test on changes in the scale.
}

\usage{
scale_cusum(x, version = "empVar", method = "kernel", control = list(), constant = 1.4826, beta = 0.5, fpc = TRUE, tol = 1e-08)
}

\arguments{
\item{x}{time series (numeric or ts vector).}

\item{version}{variance estimation method (see \code{\link{scale_stat}}).}

\item{method}{method for estimating the long run variance.}

\item{control}{a list of control parameters.}

\item{constant}{scale factor for the MAD.}

\item{beta}{quantile of the distribution function of all absolute pairwise differences used in \code{version = "QBeta"}.}

\item{fpc}{finite population correction (boolean).}

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
TO DO
}

\references{
Gerstenberger, C., Vogel, D., and Wendler, M. (2020). Tests for scale changes based on pairwise differences. Journal of the American Statistical Association, 115(531), 1336-1348.

\enc{Dürre}{Duerre}, A. (2021+). "Finite sample correction for cusum tests", \emph{unpublished manuscript}
}

\author{
Sheila \enc{Görz}{Goerz}
}

\seealso{
\code{\link{scale_stat}}, \code{\link{medianDiff}}, \code{\link{lrv}}, \code{\link{pKSdist}}
}

\examples{
TO DO
}