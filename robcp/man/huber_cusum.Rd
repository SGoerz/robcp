\name{huber_cusum}
\alias{huber_cusum}

\title{
Huberized CUSUM test
}
\description{
Performs a CUSUM test on data transformed by \code{\link{psi}}. Depending in the choosen psi-function different types of changes can be detected.
}
\usage{
huber_cusum(x, fun = "HLm", tol = 1e-8, b_n, k, constant)
}

\arguments{
  \item{x}{numeric vector containing a single time series or a numeric matrix containing multiple time series (column-wise).}
  \item{fun}{character string specifiyng the transformation function \eqn{\psi}, see details.}
  \item{tol}{tolerance of the distribution function (numeric), which is used do compute p-values.}
  \item{b_n}{bandwidth, which is used to estimate the long run variance, see the help page of \code{\link{sigma2}} for details.}
  \item{k}{numeric bound used in \code{\link{psi}}.}
	\item{constant}{scale factor of the MAD. Default is 1.4826.}
}

\value{
A list of the class "htest" containing the following compontents:
  \item{statistic}{value of the test statistic (numeric).}
  \item{p.value}{p-value (numeric).}
  \item{alternative}{alternative hypothesis (character string).}
  \item{method}{name of the performed test(character string).}
  \item{data.name}{name of the data (character string).}
}

\author{
Sheila \enc{Görz}{Goerz}
}

\details{
The function performs a Huberized CUSUM test. First the data is transformed by a suitable psi-function. To detect changes in location one can apply \code{fun = "HLm"}, \code{"HLg"}, \code{"VLm"} or \code{"VLg"}, for changes in scale \code{fun = "HCm"} is avaliable and for changes in the dependence respectively covariance structure \code{fun = "HCm"}, \code{"HCg"}, \code{"VCm"} and \code{"VCg"} are possible. The actual definitions of the psi-functions can be found in the help page of \code{\link{psi}}. Denote \eqn{Y_1,\ldots,Y_n} the transformed time series. If \eqn{Y_1} is one dimensional, then the teststatistik
\deqn{V=\max_{k=1,\ldots,n} \frac{1}{\sqrt{n}\sigma} |\sum_{i=1}^k Y_i-\frac{k}{n} \sum_{i=1}^n Y_i|
}
is calculated, where \eqn{\sigma^2} is an estimator for the long run variance, see the help function of \code{\link{sigma2}} for details. \eqn{V} is asymptoticaly Kolmogorov-Smirnov distributed. We use a finite sample correction \eqn{V+0.58/\sqrt{n}} to improve finite sample performance. \cr
If \eqn{Y[1]} is multivariate, then the test statistic
\deqn{W=\max_{k=1,\ldots,n} \frac{1}{n}(\sum_{i=1}^k Y_i-\frac{k}{n} \sum_{i=1}^n Y_i)^T\Sigma^{-1}(\sum_{i=1}^k Y_i-\frac{k}{n} \sum_{i=1}^n Y_i)
}
is computed, where \eqn{\Sigma} is the long run covariance, see also \code{\link{sigma2}} for details. \eqn{W} is asymptotically distributed like the maximum of a squared Bessel bridge. We use the identity derived in Kiefer to derive p-values. Like in the one dimensional case we use a finite sample correction \eqn{(\sqrt{W}+0.58/\sqrt{n})^2}.
}

\references{
\enc{Dürre}{Duerre}, A. and Fried, R. (2019). "Robust change point tests by bounded transformations", \url{https://arxiv.org/abs/1905.06201}

Kiefer, J. (1959). "K-sample analogues of the Kolmogorov-Smirnov and Cramer-V. Mises tests", \emph{The Annals of Mathematical Statistics}, 420--447.
}

\seealso{
\code{\link{sigma2}}, 
\code{\link{psi}}, \code{\link{h_cumsum}}, \code{\link{teststat}}, 
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
