\name{cor_stat}
\alias{cor_stat}

\title{Test statistic to detect Fluctuation Changes}

\description{
Computes the test statistic for a CUSUM-based tests on changes in Spearman's rho.
}

\usage{
cor_stat(x, version = "rho", control = list())
}

\arguments{
\item{x}{time series (numeric or ts vector).}
\item{version}{version of the test. Either \code{rho} or \code{tau}.}
\item{control}{a list of control parameters.}
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
Wied, D., Dehling, H., Van Kampen, M., and Vogel, D. (2014). A fluctuation test for constant Spearman’s rho with nuisance-free limit distribution. Computational Statistics & Data Analysis, 76, 723-736.
}

\seealso{
\code{\link{lrv}}, \code{\link{cor_cusum}}
}

\examples{
TO DO
}
