\encoding{utf8}
\name{pKSdist}
\alias{pKSdist}
\alias{pBessel}

\title{
Asymptotic cumulative distribution for the CUSUM Test statistic
}
\description{
Computes the asymptotic cumulative distribution of the statistic of \code{\link{CUSUM}}.
}
\usage{
pKSdist(tn, tol = 1e-8)
pBessel(tn, p)
}

\arguments{
  \item{tn}{vector of test statistics (numeric). For \code{pBessel} length of \code{tn} has to be 1.}
  \item{p}{dimension of time series (integer). If \code{p} is equal to 1 \code{pBessel} uses \code{pKSdist} to compute the corresponding probability.}
  \item{tol}{tolerance (numeric).}
}

\details{
For a single time series, the distribution is the same distribution as in the two sample Kolmogorov Smirnov Test, namely the distribution of the maximal value of the absolute values of a Brownian bridge. It is computated as follows (Durbin, 1973 and van Mulbregt, 2018):

For \eqn{t_n(x) < 1}:

\deqn{
P(t_n(X) \le t_n(x)) = 
\frac{\sqrt{2 \pi}}{t_n(x)} t (1 + t^8(1 + t^{16}(1 + t^{24}(1 + ...))))}
up to \eqn{t^{8 k_{max}},
k_{max} = \lfloor \sqrt{2 - \log(tol)}\rfloor}, where \eqn{t = \exp(-\pi^2 / (8x^2))}

else:

\deqn{
P(t_n(X) \le t_n(x)) = 2  \sum_{k = 1}^{\infty} (-1)^{k - 1} \exp(-2 k^2 x^2)}
until
\eqn{|2 (-1)^{k - 1} \exp(-2 k^2 x^2) - 2 (-1)^{(k-1) - 1}  \exp(-2 (k-1)^2 x^2)| \le tol.
}

In case of multiple time series, the distribution equals that of the maximum of an \code{p} dimensional squared Bessel bridge. It can be computed by

\deqn{P(t_n(X) \le t_n(x)) = 
\frac{4}{ \Gamma(p / 2) 2^{p / 2} t_n^p }  \sum_{i = 1}^{\infty}  \frac{(\gamma_{(p - 2)/2, n})^{p - 2} \exp(-(\gamma_{(p - 2)/2, n})^2 / (2t_n^2))}{J_{p/2}(\gamma_{(p - 2)/2, n})^2 }}
(Kiefer, 1959), where \eqn{J_p} is the Bessel function of first kind and p-th order, \eqn{\Gamma} is the gamma function and \eqn{\gamma_{p, n}} denotes the n-th zero of \eqn{J_p}.
}

\value{
vector of \eqn{P(t_n(X) \le tn[i])}.
}

\references{
Durbin, James. (1973) "Distribution theory for tests based on the sample distribution function." \emph{Society for Industrial and Applied Mathematics}.

van Mulbregt, P. (2018) "Computing the Cumulative Distribution Function and Quantiles of the limit of the Two-sided Kolmogorov-Smirnov Statistic." arXiv preprint arXiv:1803.00426.

/src/library/stats/src/ks.c rev60573

Kiefer, J. (1959). "K-sample analogues of the Kolmogorov-Smirnov and Cramer-V. Mises tests", \emph{The Annals of Mathematical Statistics}, 420--447.
}

\author{
Sheila \enc{Görz}{Goerz}, Alexander \enc{Dürre}{Duerre}
}

\seealso{
\code{\link{psi}}, \code{\link{CUSUM}}, \code{\link{psi_cumsum}}, \code{\link{huber_cusum}}
}

\examples{
# single time series
timeSeries <- c(rnorm(20, 0), rnorm(20, 2))
tn <- CUSUM(timeSeries)

pKSdist(tn)

# two time series
timeSeries <- matrix(c(rnorm(20, 0), rnorm(20, 2), rnorm(20, 1), rnorm(20, 3)), 
                     ncol = 2)
tn <- CUSUM(timeSeries)

pBessel(tn, 2)

}

