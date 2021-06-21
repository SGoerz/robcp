\name{kthPair}
\alias{kthPair}

\title{
K-th largest element in a sum of sets.
}
\description{
Selects the k-th largest element of X + Y, a sum of sets. X + Y denotes the set \eqn{\{x + y | x \in X, y \in Y\}}. 
}
\usage{
kthPair(x, y, k, k2 = NA)
}

\arguments{
  \item{x}{Numeric vector.}
  \item{y}{Numeric vector.}
  \item{k}{Index of element to be selected. Must be an integer and between 1 and the product of the lengths of x and y.}
  \item{k2}{Optional second index. \code{k} and \code{k2} must be consecutive. Useful, if X + Y has an even number of elements and the median has to be calculated.}
}

\details{
...missing...
}

\value{
K-th largest value (numeric).
}

\references{
Johnson, D. B., & Mizoguchi, T. (1978). Selecting the K-th Element in X+Y and X_1+X_2+ ... +X_m. SIAM Journal on Computing, 7(2), 147-153.
}

\author{
Sheila \enc{Görz}{Goerz}
}

\examples{
set.seed(1895)
x <- rnorm(100)
y <- runif(100)

kthPair(x, y, 5000)
}
