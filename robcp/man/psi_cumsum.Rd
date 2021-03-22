\encoding{utf8}
\name{psi_cumsum}
\alias{psi_cumsum}

\title{
Cumulative sum of transformed vectors
}
\description{
Computes the cumulative sum of a transformed numeric vector or matrix. Transformation function is \code{\link{psi}}.
}

\usage{
psi_cumsum(y, fun = "HLm", k, constant)
}

\arguments{
	\item{y}{numeric vector containing a single time series or a numeric matrix containing multiple time series (column-wise).}
	\item{fun}{character string specifiyng the transformation function \eqn{\psi}.}
	\item{k}{numeric bound used in \code{\link{psi}}.}
	\item{constant}{scale factor of the MAD. Default is 1.4826.}
}

\value{
Numeric vector or matrix containing the cumulative sums of the transformed values. In case of a matrix, cumulative sums are being computed for each time series (column) independentely.
}

\details{
Prior to computing the sums, y is being transformed by the function \code{fun}.
}

\seealso{
\code{\link{psi}}.
}

\author{
Sheila \enc{Görz}{Goerz}
}


\examples{
psi_cumsum(rnorm(100))
}
