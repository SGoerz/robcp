\name{scale_stat}
\alias{scale_stat}

\title{Test statistic to detect Scale Changes}

\description{
Computes the test statistic for CUSUM-based tests on scale changes.
}

\usage{
scale_stat(x, version = "empVar", method = "kernel", control = list(), constant = 1.4826, beta = 0.5)
}

\arguments{
\item{x}{time series (numeric or ts vector).}

\item{version}{variance estimation method. One of \code{"empVar"}, \code{"MD"}, \code{"GMD"}, \code{"MAD"}, \code{"QBeta"}.}

\item{method}{method for estimating the long run variance.}

\item{control}{a list of control parameters.}

\item{constant}{scale factor for the MAD.}

\item{beta}{quantile of the distribution function of all absolute pairwise differences used in \code{version = "QBeta"}.}
}

\value{
Test statistic (numeric value) with the attribute cp-location 
       indicating at which index a change point is most likely. Is an S3 
       object of the class cpStat
}

\details{
TO DO
}

\author{
Sheila \enc{Görz}{Goerz}
}
\references{
Gerstenberger, C., Vogel, D., and Wendler, M. (2020). Tests for scale changes based on pairwise differences. Journal of the American Statistical Association, 115(531), 1336-1348.
}

\seealso{
\code{\link{medianDiff}}, \code{\link{lrv}}
}

\examples{
TO DO
}
