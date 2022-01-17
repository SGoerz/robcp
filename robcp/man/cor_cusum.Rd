\name{cor_cusum}
\alias{cor_cusum}
\title{A fluctuation test for constant Spearman’s rho with nuisance-free limit distribution}

\description{
Performs a CUSUM-based test on changes in Spearman's rho.
}

\usage{
cor_cusum(x, version = "rho", control = list(), fpc = TRUE, tol = 1e-08)
}

\arguments{
\item{x}{time series (numeric or ts vector).}
\item{version}{version of the test. Either \code{rho} or \code{tau}.}
\item{control}{a list of control parameters.}
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
Wied, D., Dehling, H., Van Kampen, M., and Vogel, D. (2014). A fluctuation test for constant Spearman’s rho with nuisance-free limit distribution. Computational Statistics & Data Analysis, 76, 723-736.

\enc{Dürre}{Duerre}, A. (2021+). "Finite sample correction for cusum tests", \emph{unpublished manuscript}
}

\author{
Sheila \enc{Görz}{Goerz}
}

\seealso{
\code{\link{cor_stat}}, \code{\link{lrv}}, \code{\link{pKSdist}}
}

\examples{
TO DO
}
