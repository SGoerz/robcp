\name{huber_cusum}
\alias{huber_cusum}

\title{
Huberized CUSUM test
}
\description{
Performs a CUSUM test on data transformed by \code{\link{psi}}. Depending on the chosen psi-function different types of changes can be detected.
}
\usage{
huber_cusum(x, fun = "HLm", k, constant = 1.4826, method = "kernel",
            control = list(), fpc = TRUE, tol = 1e-8, ...)
}

\arguments{
  \item{x}{numeric vector containing a single time series or a numeric matrix containing multiple time series (column-wise).}
  \item{fun}{character string specifying the transformation function \eqn{\psi}, see details.}
  \item{k}{numeric bound used in \code{\link{psi}}.}
	\item{constant}{scale factor of the MAD. Default is 1.4826.}
	\item{method}{method for estimating the long run variance.}
	\item{control}{a list of control parameters for the estimation of the long run variance (cf. \code{\link{lrv}}).}
	\item{fpc}{finite population correction (boolean).}
	\item{tol}{tolerance of the distribution function (numeric), which is used to compute p-values.}
	\item{...}{further arguments to be passed to \code{\link{CUSUM}}.}
}

\value{
A list of the class "htest" containing the following components:
  \item{statistic}{value of the test statistic (numeric).}
  \item{p.value}{p-value (numeric).}
  \item{alternative}{alternative hypothesis (character string).}
  \item{method}{name of the performed test (character string).}
  \item{cp.location}{index of the estimated change point location (integer).}
  \item{data.name}{name of the data (character string).}
}

\author{
Sheila \enc{Görz}{Goerz}
}

\details{
The function performs a Huberized CUSUM test. It tests the null hypothesis \eqn{H_0: \boldsymbol{\theta}} does not change for \eqn{x} against the alternative of a change, where \eqn{\boldsymbol{\theta}} is the parameter vector of interest. \eqn{k} is called a 'change point'. First the data is transformed by a suitable psi-function. To detect changes in location one can apply \code{fun = "HLm"}, \code{"HLg"}, \code{"SLm"} or \code{"SLg"} and the hypothesis pair is
\deqn{H_0: \mu_1 = ... = \mu_n}
\deqn{vs.}
\deqn{H_1: \exists k \in \{1, ..., n-1\}: \mu_k \neq \mu_{k+1}}
where \eqn{\mu_t = E(X_t)} and \eqn{n} is the length of the time series. For changes in scale \code{fun = "HCm"} is available and for changes in the dependence respectively covariance structure \code{fun = "HCm"}, \code{"HCg"}, \code{"SCm"} and \code{"SCg"} are possible. The hypothesis pair is the same as in the location case, only with \eqn{\mu_i} being replaced by \eqn{\Sigma_i}, \eqn{\Sigma_i = Cov(X_i)}. Exact definitions of the psi-functions can be found on the help page of \code{\link{psi}}. 

Denote by \eqn{Y_1,\ldots,Y_n} the transformed time series. If \eqn{Y_1} is one-dimensional, then the test statistic
\deqn{V_n = \max_{k=1,\ldots,n} \frac{1}{\sqrt{n}\sigma} \left|\sum_{i=1}^k Y_i-\frac{k}{n} \sum_{i=1}^n Y_i\right|
}
is calculated, where \eqn{\sigma^2} is an estimator for the long run variance, see the help function of \code{\link{lrv}} for details. \eqn{V} is asymptotically Kolmogorov-Smirnov distributed. If \code{fpc} is \code{TRUE} we use a finite population correction \eqn{V+0.58/\sqrt{n}} to improve finite sample performance.  \cr
If \eqn{Y_1} is multivariate, then the test statistic
\deqn{W_n=\max_{k=1,\ldots,n} \frac{1}{n}\left(\sum_{i=1}^k Y_i-\frac{k}{n} \sum_{i=1}^n Y_i\right)' \Sigma^{-1}\left(\sum_{i=1}^k Y_i-\frac{k}{n} \sum_{i=1}^n Y_i\right)
}
is computed, where \eqn{\Sigma} is the long run covariance, see also \code{\link{lrv}} for details. \eqn{W} is asymptotically distributed like the maximum of a squared Bessel bridge. We use the identity derived by Kiefer (1959) to derive p-values. Like in the one dimensional case if \code{fpc} is \code{TRUE} we use a finite sample correction \eqn{(\sqrt{W}+0.58/\sqrt{n})^2}.

The change point location is estimated as the time point \eqn{k} for which the CUSUM process takes its maximum. 
}

\references{
\enc{Dürre}{Duerre}, A. and Fried, R. (2019). "Robust change point tests by bounded transformations", \url{https://arxiv.org/abs/1905.06201}

Kiefer, J. (1959). "K-sample analogues of the Kolmogorov-Smirnov and Cramer-V. Mises tests", \emph{The Annals of Mathematical Statistics}, 420--447.

}

\seealso{
\code{\link{lrv}}, 
\code{\link{psi}}, \code{\link{psi_cumsum}}, \code{\link{CUSUM}}, 
\code{\link{pKSdist}}
}

\examples{
set.seed(1895)

#time series with a structural break at t = 20
Z <- c(rnorm(20, 0), rnorm(20, 2))
huber_cusum(Z) 

# two time series with a structural break at t = 20
timeSeries <- matrix(c(rnorm(20, 0), rnorm(20, 2), rnorm(20, 1), rnorm(20, 3), 
                     ncol = 2))
                     
huber_cusum(timeSeries)
}
