\encoding{utf8}
\name{teststat}
\alias{teststat}

\title{
Test statistic for the Huberized CUSUM Test
}

\description{
Computes the test statistic for a structural break test called 'Huberized CUSUM Test'.
}

\usage{
teststat(y, fun = "HLm", b_n, k, constant)
}

\arguments{
	\item{y}{vector or matrix with each column representing a time series (numeric).}
	\item{fun}{character string specifiyng the transformation function \eqn{\psi}.}
	\item{b_n}{for \code{sigma2}.}
	\item{k}{numeric bound used in \code{\link{psi}}.}
	\item{constant}{scale factor of the MAD. Default is 1.4826.}
}

\value{
test statistic (numeric value).
}

\details{
y is transformed by \code{fun}. Let x be the resulting vector or matrix and n be the length of a time series. 

In case of a vector the test statistic can be written as \deqn{max_{k = 1, ..., n}\frac{1}{\sqrt{n} \sigma}|\sum_{i = 1}^{k} x_i - (k / n) \sum_{i = 1}^n x_i|,} where \eqn{\sigma} is the square root of \code{sigma2}.

In case of a matrix the test statistic follows as
\deqn{max_{k = 1, ..., n}\frac{1}{n}(\sum_{i = 1}^{k} X_i - \frac{k}{n} \sum_{i = 1}^{n} X_i)^T \Sigma^{-1} (\sum_{i = 1}^{k} X_i - \frac{k}{n} \sum_{i = 1}^{n} X_i),} where \eqn{X_i} denotes the ith row of x and \eqn{\Sigma^{-1}} is the inverse of \code{sigma2}.
}

\seealso{
\code{\link{h_cumsum}}, 
\code{\link{psi}}
}


\author{
Sheila \enc{Görz}{Goerz}
}


\examples{
# time series with structural break at t = 20
ts <- c(rnorm(20, 0), rnorm(20, 2))

teststat(ts)
}
