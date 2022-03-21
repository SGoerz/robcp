\name{Qalpha}
\alias{Qalpha}

\title{\eqn{Q^{\alpha}}}

\description{
Estimates Q-alpha using the first \eqn{k} elements of a time series.
}

\usage{
Qalpha(x, alpha)
}

\arguments{
\item{x}{numeric vector.}
\item{alpha}{quantile. Numeric value in (0, 1]}
}

\details{
\deqn{\hat{Q}^{\alpha}_{1:k} = U_{1:k}^{-1}(\alpha) = \inf\{x | \alpha \leq U_{1:k}(x)\},}
where \eqn{U_{1:k}} is the empirical distribtion function of \eqn{|x_i - x_j|, \, 1 \leq i < j \leq k}. 
}

\value{
Numeric vector of \eqn{Q^{\alpha}}-s estimated using \code{x[1], ..., x[k]}, \eqn{k = 1, ..., n}, \eqn{n} being the length of \code{x}.
}

\examples{
x <- rnorm(10)
Qalpha(x, 0.5)
}
