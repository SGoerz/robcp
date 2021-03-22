\name{modifChol}
\alias{modifChol}

\title{
Revised Modified Cholesky Factorization
}

\description{
Computes the revised modified Cholesky factorization described in Schnabel and Eskow (1999).
}
\usage{
modifChol(x, tau = .Machine$double.eps^(1 / 3), 
           tau_bar = .Machine$double.eps^(2 / 3), mu = 0.1)
}

\arguments{
  \item{x}{a symmetric matrix.}
  \item{tau}{(machine epsilon)^(1/3).}
  \item{tau_bar}{(machine epsilon^(2/3)).}
  \item{mu}{numeric, \eqn{0 < \mu \le 1}.}
}

\details{
\code{modif.chol} computes the revised modified Cholesky Factorization of a symmetric, not neccessarily positive definite matrix x + E such that \eqn{L'L = x + E} for \eqn{E \ge 0}.
}

\value{
Upper triangular matrix \eqn{L} of the form \eqn{L'L = x + E}.
The attribute \code{swaps} is a vector of the lenght of dimension of x. It cointains the indices of the rows and columns that were swapped in x in order to compute the modified Cholesky factorization. For example if the \code{i}-th element of \code{swaps} is the number \code{j}, then the \code{i}-th and the \code{j}-th row and column were swapped. To reconstruct the original matrix swaps has to be read backwards.
}

\references{
Schnabel, R. B., & Eskow, E. (1999). "A revised modified Cholesky factorization algorithm" SIAM Journal on optimization, 9(4), 1135-1148.
}

\author{
Sheila Görz
}

\examples{
y <- matrix(runif(9), ncol = 3)
x <- psi(y)
modifChol(lrv(x))
}
