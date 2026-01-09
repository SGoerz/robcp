\name{cor_cusum}
\alias{cor_cusum}
\title{A CUSUM-type test to detect changes in the correlation.}

\description{
Performs a CUSUM-based test on changes in Spearman's rho or Kendall's tau.
}

\usage{
cor_cusum(x, version = c("tau", "rho"), method = "kernel", control = list(), 
          fpc = TRUE, tol = 1e-08, plot = FALSE)
}

\arguments{
\item{x}{time series (matrix or ts object with numeric/integer values).}
\item{version}{version of the test. Either \code{rho} or \code{tau}.}
\item{method}{method for estimating the long run variance.}
\item{control}{a list of control parameters.}
\item{fpc}{finite population correction (boolean).}
\item{tol}{tolerance of the distribution function (numeric), which is used do compute p-values.}
\item{plot}{should the test statistic be plotted (cf. \code{\link{plot.cpStat}})? Boolean.}
}

\value{
A list of the class "htest" containing the following components:
  \item{statistic}{return value of the function \code{\link{cor_stat}}.}
  \item{p.value}{p-value (numeric).}
  \item{alternative}{alternative hypothesis (character string).}
  \item{method}{name of the performed test (character string).}
  \item{cp.location}{index of the estimated change point location (integer).}
  \item{data.name}{name of the data (character string).}
  \item{lrv}{list containing the compontents \code{method}, \code{param} and \code{value}.}
}


\details{
The function perform a CUSUM-type test on changes in the correlation of a time series \eqn{x}. Formally, the hypothesis pair can be written as
\deqn{H_0: \xi_1 = ... = \xi_n}
\deqn{vs.}
\deqn{H_1: \exists k \in \{1, ..., n-1\}: \xi_k \neq \xi_{k+1}}
where \eqn{\xi_i} is a fluctuation measure (either Spearman's rho or Kendall's tau) for the \eqn{i}-th observations and \eqn{n} is the length of the time series. \eqn{k} is called a 'change point'.

The test statistic is computed using \code{\link{cor_stat}} and asymptotically follows a Kolmogorov distribution. To derive the p-value, the funtion \code{\link{pKSdist}} is used.
}

\references{
Wied, D., Dehling, H., Van Kampen, M., and Vogel, D. (2014). A fluctuation test for constant Spearman’s rho with nuisance-free limit distribution. Computational Statistics & Data Analysis, 76, 723-736.

\enc{Dürre}{Duerre}, A. (2022+). "Finite sample correction for cusum tests", \emph{unpublished manuscript}
}

\author{
Sheila \enc{Görz}{Goerz}
}

\seealso{
\code{\link{cor_stat}}, \code{\link{lrv}}, \code{\link{pKSdist}}
}

\examples{
### first: generate a time series with a burn-in period of m and a change point k = n/2
require(mvtnorm)
n <- 500
m <- 100
N <- n + m
k <- m + floor(n * 0.5)
n1 <- N - k

## Spearman's rho:
rho <- c(0.4, -0.9)
  
# serial dependence:
theta1 <- 0.3
theta2 <- 0.2
theta <- cbind(c(theta1, 0), c(0, theta2))
q <- rho * sqrt( (theta1^2 + 1) * (theta2^2 + 1) / (theta1 * theta2 + 1))
# shape matrices of the innovations:
S0 <- cbind(c(1, q[1]), c(q[1], 1))
S1 <- cbind(c(1, q[2]), c(q[2], 1))

e0 <- rmvt(k, S0, 5)
e1 <- rmvt(n1, S1, 5)
e <- rbind(e0, e1)
# generate the data:
x <- matrix(numeric(N * 2), ncol = 2)
x[1, ] <- e[1, ]
invisible(sapply(2:N, function(i) x[i, ] <<- e[i, ] + theta \%*\% e[i-1, ]))
x <- x[-(1:m), ]

cor_cusum(x, "rho")


## Kendall's tau
S0 <- cbind(c(1, rho[1]), c(rho[1], 1))
S1 <- cbind(c(1, rho[2]), c(rho[2], 1))
e0 <- rmvt(k, S0, 5)
e1 <- rmvt(n1, S1, 5)
e <- rbind(e0, e1)
x <- matrix(numeric(N * 2), ncol = 2)
x[1, ] <- e[1, ]
# AR(1):
invisible(sapply(2:N, function(i) x[i, ] <<- 0.8 * x[i-1, ] + e[i, ]))
x <- x[-(1:m), ]

cor_cusum(x, version = "tau")

}
