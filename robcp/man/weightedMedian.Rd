\name{weightedMedian}
\alias{weightedMedian}
\title{Weighted Median}
\description{Computes the weighted median of a numeric vector.}
\usage{
weightedMedian(x, w)
}
\arguments{
\item{x}{Numeric vector.}

\item{w}{Integer vector of weights.}
}
\value{
Weighted median of x with respect to w.
}
\details{
Here, the median of an even length \eqn{n} of \eqn{x} is defined as \eqn{x_{(n/2 + 1)}} if \eqn{x_{(i)}} is the \eqn{i}-th largest element in \eqn{x}.
}

\examples{
x <- c(1, 4, 9)
w <- c(5, 1, 1)
weightedMedian(x, w)
}
