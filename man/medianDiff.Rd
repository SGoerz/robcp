\name{medianDiff}
\alias{medianDiff}
\title{Median of the set X - Y}
\usage{
medianDiff(x, y)
}
\arguments{
\item{x, y}{Numeric vectors}
}
\value{
The median element of X - Y
}
\description{
Computes the median of the set X - Y. X - Y denotes the set \eqn{\{x - y | x \in X, y \in Y}\}.
}

\details{
Special case of the function \code{\link{kthPair}}.
}

\author{
Sheila \enc{Görz}{Goerz}
}
\references{
Johnson, D. B., & Mizoguchi, T. (1978). Selecting the K-th Element in X+Y and X_1+X_2+ ... +X_m. SIAM Journal on Computing, 7(2), 147-153.
}
\examples{
x <- rnorm(100)
y <- runif(100)
medianDiff(x, y)
}
