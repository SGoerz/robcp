\name{lrv}
\alias{lrv}

\title{
Long Run Variance
}
\description{
Estimates the long run variance respectively covariance matrix of the supplied time series.
}

\usage{
lrv(x, method = c("kernel", "subsampling", "none"), control = list())
}

\arguments{
  \item{x}{vector or matrix with each column representing a time series (numeric).}
  \item{method}{method of estimation. Options are \code{kernel}, \code{subsampling}, \code{bootstrap} and \code{none}.}
  \item{control}{a list of control parameters. See 'Details'.}
}

\details{
The long run variance equals the limit of \eqn{n} times the variance of the arithmetic mean of a short range dependent time series, where \eqn{n} is the length of the time series. It is used to standardize tests concering the mean on dependent data.

If \code{method = "none"}, no long run variance estimation is performed and the value 1 is returned (i.e. it does not alterate the test statistic). 

The \code{control} argument is a list that can supply any of the following components:
\describe{
\item{\code{kFun}}{Kernel function (character string). More in 'Notes'.}
\item{\code{b_n}}{Bandwidth (numeric > 0 and smaller than sample size).}
\item{\code{gamma0}}{Only use estimated variance if estimated long run variance is < 0? Boolean.}
\item{\code{l}}{Block length (numeric > 0 and smaller than sample size).}
\item{\code{overlapping}}{Overlapping subsampling estimation? Boolean.}
\item{\code{distr}}{Tranform observations by their empirical distribution function? Boolean. Default is \code{FALSE}.}
\item{\code{B}}{Bootstrap repetitions (integer).}
\item{\code{seed}}{RNG seed (numeric).}
\item{\code{version}}{What property does the CUSUM test test for? Character string, details below.}
\item{\code{loc}}{Estimated location corresponding to \code{version}. Numeric value, details below.}
\item{\code{scale}}{Estimated scale corresponding to \code{version}. Numeric value, details below.}
}

\strong{Kernel-based estimation}

The kernel-based long run variance estimation is available for various testing scenarios (set by \code{control$version}) and both for one- and multi-dimensional data. It uses the bandwidth \eqn{b_n = } \code{control$b_n} and kernel function \eqn{k(x) = } \code{control$kFun}. For tests on certain properties also a corresponding location \code{control$loc} (\eqn{m_n}) and scale \code{control$scale} (\eqn{v_n}) estimation needs to be supplied. Supported testing scenarios are:

\itemize{
  \item \code{"mean"}
  \itemize{
    \item 1-dim. data:
    \deqn{\hat{\sigma}^2 = \frac{1}{n} \sum_{i = 1}^n (x_i - \bar{x})^2 + \frac{2}{n} \sum_{h = 1}^{b_n} \sum_{i = 1}^{n - h} (x_i - \bar{x}) (x_{i + h} - \bar{x}) k(h / b_n).}
    If \code{control$distr = TRUE}, then the long run variance is estimated on the empirical distribution of \eqn{x}. The resulting value is then multiplied with \eqn{\sqrt{\pi} / 2}.
    
    Default values: \code{b_n} = \eqn{0.9 n^{1/3}}, \code{kFun = "bartlett"}.
    
    \item multivariate time series:
    The \eqn{k,l}-element of \eqn{\Sigma} is estimated by
    \deqn{\hat{\Sigma}^{(k,l)} = \frac{1}{n} \sum_{i,j = 1}^{n}(x_i^{(k)} - \bar{x}^{(k)})  (x_j^{(l)} -        \bar{x}^{(l)})  k((i-j) / b_n),}
    \eqn{k, l = 1, ..., m}.

    Default values: \code{b_n} = \eqn{\log_{1.8 + m / 40}(n / 50)}, \code{kFun = "bartlett"}.
  }
  
  \item \code{"empVar"} for tests on changes in the empirical variance.
  \deqn{\hat{\sigma}^2 = \sum_{h = -(n-1)}^{n-1} W \left( \frac{|h|}{b_n} \right) \frac{1}{n} \sum_{i = 1}^{n - |h|} ((x_i - m_n)^2 - v_n)((x_{i+|h|} - m_n)^2 - v_n).}
  
  Default values: \eqn{m_n =} \code{mean(x)}, \eqn{v_n = } \code{var(x)}.
  
  \item \code{"MD"} for tests on a change in the median deviation.
  \deqn{\hat{\sigma}^2 = \sum_{h = -(n-1)}^{n-1} W \left( \frac{|h|}{b_n} \right) \frac{1}{n} \sum_{i = 1}^{n - |h|} (|x_i - m_n| - v_n)(|x_{i+|h|} - m_n| - v_n).}
  
  Default values: \eqn{m_n =} \code{median(x)}, \eqn{v_n = \frac{1}{n-1} \sum_{i = 1}^n |x_i - m_n|}.
  
  \item \code{"GMD"} for tests on changes in Gini's mean difference.
  \deqn{\hat{\sigma}^2 = 4 \sum_{h = -(n-1)}^{n-1} W \left( \frac{|h|}{b_n} \right) \frac{1}{n} \sum_{i = 1}^{n - |h|} \hat{\phi}_n(x_i)\hat{\phi}_n(x_{i+|h|})}
  with \eqn{\hat{\phi}_n(x) = n^{-1} \sum_{i = 1}^n |x - x_i| - v_n}.
  
  Default value: \eqn{v_n =} \eqn{\frac{2}{n(n-1)} \sum_{1 \leq i < j \leq n} |x_i - x_j|.}
  
  %\item \code{"MAD"} for tests on changes in the median absolute deviation (from the median).
  %\deqn{\hat{\sigma}^2 = \hat{f}_Z(v_n)^{-1} \sum_{h = -(n-1)}^{n-1} W \left( \frac{|h|}{b_n} \right) \frac{1}{n} \sum_{i = 1}^{n - |h|} \hat{\phi}_n(x_i)\hat{\phi}_n(x_{i+|h|}),}
   %where \eqn{\hat{\phi}_n(x) = 1_{\{|x - m_n| \leq v_n\}} - 0.5} and 
   %\deqn{\hat{f}_Z(t) = \frac{1}{n h_n} \sum_{i = 1}^n K\left(\frac{|x_i - m_n| - t}{h_n}\right)} 
   %the kernel density estimation of the density of \eqn{|X - med(F^X)|}, \eqn{h_n =} \code{IQR(x)*n^{-1/3}}, and \eqn{K} is the quatratic kernel function.
   %
   %Default values: \eqn{m_n =} \code{median(x)}, \eqn{v_n =} \code{mad(x)}.
  
  \item \code{"Qalpha"} for tests on changes in \code{\link{Qalpha}}.
  \deqn{\hat{\sigma}^2 = \frac{4}{\hat{u}(v_n)} \sum_{h = -(n-1)}^{n-1} W \left( \frac{|h|}{b_n} \right) \frac{1}{n} \sum_{i = 1}^{n - |h|} \hat{\phi}_n(x_i)\hat{\phi}_n(x_{i+|h|}),}
  where \eqn{\hat{\phi}_n(x) = n^{-1} \sum_{i = 1}^n 1_{\{|x - x_i| \leq v_n\}} - m_n} and 
  \deqn{\hat{u}(t) = \frac{2}{n(n-1)h_n} \sum_{1 \leq i < j \leq n} K\left(\frac{|x_i - x_j| - t}{h_n}\right)}
  the kernel density estimation of the densitiy \eqn{u} corresponding to the distribution function \eqn{U(t) = P(|X-Y| \leq t)}, \eqn{h_n =} \code{IQR(x)}\eqn{n^{-\frac{1}{3}}} and \eqn{K} is the quatratic kernel function.
  
  Default values: \eqn{m_n = \alpha = 0.5}, \eqn{v_n =} \code{Qalpha(x, m_n)[n-1]}.
  
  \item \code{"tau"} for tests in changes in Kendall's tau.
  
  Only available for bivariate data: assume that the given data \code{x} has the format \eqn{(x_i, y_i)_{i = 1, ..., n}}.
  \deqn{\hat{\sigma}^2 = \sum_{h = -(n-1)}^{n-1} W \left( \frac{|h|}{b_n} \right) \frac{1}{n} \sum_{i = 1}^{n - |h|} \hat{\phi}_n((x_i, y_i))\hat{\phi}_n((x_{i+|h|}, y_{i+|h|}),}
  where \eqn{\hat{\phi}_n(x) = 4 F_n(x, y) - 2F_{X,n}(x) 2 - F_{Y,n}(y) + 1 - v_n} and \eqn{F_n}, \eqn{F_{X,n}} and \eqn{F_{Y,n}} are the empirical distribution functions of \eqn{((X_i, Y_i))_{i = 1, ..., n}}, \eqn{(X_i)_{i = 1, ..., n}} and \eqn{(Y_i)_{i = 1, ..., n}}.
  
  Default value: \eqn{v_n = \hat{\tau}_n = \frac{2}{n(n-1)} \sum_{1 \leq i < j \leq n} sign\left((x_j - x_i)(y_j - y_i)\right)}.
  
  \item \code{"rho"} for tests on changes in Spearman's rho.
  
  Only availabe for \eqn{d}-variate data with \eqn{d > 1}: assume that the given data \code{x} has the format \eqn{(x_{i,j} | i = 1, ..., n; j = 1, ..., d)}.
  \deqn{\hat{\sigma}^2 = a(d)^2 2^{2d} \left\{ \sum_{h = -(n-1)}^{n-1} K\left( \frac{|h|}{b_n} \right) \left( \sum_{i = 1}^{n-|h|} n^{-1} \prod_{j = 1}^d \hat{\phi}_n(x_i, x_j) \hat{\phi}_n(x_{i+|h|}, x_j) - M^2 \right) \right\} ,}
   where \eqn{a(d) = (d+1) / (2^d - d - 1)},  \eqn{M = n^{-1} \sum_{i = 1}^n \prod_{j = 1}^d \hat{\phi}_n(x_i, x_j)} and \eqn{\hat{\phi}_n(x, y) = 1 - \hat{U}_n(x, y)}, \eqn{\hat{U}_n(x, y) = n^{-1}} (rank of \eqn{x_{i,j}} in \eqn{x_{i,1}, ..., x_{i,n})}.
}

When \code{control$gamma0 = TRUE} (default) then negative estimates of the long run variance are replaced by the autocovariance at lag 0 (= ordinary variance of the data). The function will then throw a warning.


\strong{Subsampling estimation}

For \code{method = "subsampling"} there are an overlapping and a non-overlapping version (parameter \code{control$overlapping}). Also it can be specified if the observations x were transformed by their empirical distribution function \eqn{\tilde{F}_n} (parameter \code{control$distr}). Via \code{control$l} the block length \eqn{l} can be controlled.

If \code{control$overlapping = TRUE} and \code{control$distr = TRUE}:
\deqn{\hat{\sigma}_n = \frac{\sqrt{\pi}}{\sqrt{2l}(n - l + 1)} \sum_{i = 0}^{n-l} \left| \sum_{j = i+1}^{i+l} (F_n(x_j) - 0.5) \right|.}

Otherwise, if \code{control$distr = FALSE}, the estimator is
\deqn{\hat{\sigma}^2 = \frac{1}{l (n - l + 1)} \sum_{i = 0}^{n-l} \left( \sum_{j = i + 1}^{i+l} x_j - \frac{l}{n} \sum_{j = 1}^n x_j \right)^2.}

If \code{control$overlapping = FALSE} and \code{control$distr = TRUE}:
\deqn{\hat{\sigma} = \frac{1}{n/l} \sqrt{\pi/2} \sum_{i = 1}{n/l} \frac{1}{\sqrt{l}} \left| \sum_{j = (i-1)l + 1}^{il} F_n(x_j) - \frac{l}{n} \sum_{j = 1}^n F_n(x_j) \right|.}

Otherwise, if \code{control$distr = FALSE}, the estimator is
\deqn{\hat{\sigma}^2 = \frac{1}{n/l} \sum_{i = 1}^{n/l} \frac{1}{l} \left(\sum_{j = (i-1)l + 1}^{il} x_j - \frac{l}{n} \sum_{j = 1}^n x_j\right)^2.}


Default values: overlapping = TRUE, the block length is chosen adaptively: 
\deqn{l_n = \max{\left\{ \left\lceil n^{1/3} \left( \frac{2 \rho}{1 - \rho^2} \right)^{(2/3)} \right\rceil, 1 \right\}}}
where \eqn{\rho} is the Spearman autocorrelation at lag 1.


%\strong{Bootstrap estimation}

%If \code{method = "bootstrap"} a dependent wild bootstrap with the parameters \eqn{B = } \code{control$B}, \eqn{l = } \code{control$l} and \eqn{k(x) = } \code{control$kFun} is performed:
%\deqn{ \hat{\sigma}^2 = \sqrt{n} Var(\bar{x^*_k} - \bar{x}), k = 1, ..., B}
%A single \eqn{x_{ik}^*} is generated by \eqn{x_i^* = \bar{x} + (x_i - \bar{x}) a_i} where \eqn{a_i} are independent from the data \code{x} and are generated from a multivariate normal distribution with \eqn{E(A_i) = 0}, \eqn{Var(A_i) = 1} and \eqn{Cov(A_i, A_j) = k\left(\frac{i - j}{l}\right), i = 1, ..., n; j \neq i}. Via \code{control$seed} a seed can optionally be specified (cf. \code{\link{set.seed}}). Only \code{"bartlett"}, \code{"parzen"} and \code{"QS"} are supported as kernel functions. Uses the function \code{\link[pracma]{sqrtm}} from package \code{pracma}.

%Default values: \code{B} = 1000, \code{kFun = "bartlett"}, \code{l} is the same as for subsampling.
}


\note{
Kernel functions

\code{bartlett}:
\deqn{k(x) = (1 - |x|) * 1_{\{|x| < 1\}}}

\code{FT}:
\deqn{k(x) = 1 * 1_{\{|x| \leq 0.5\}} + (2 - 2 * |x|) * 1_{\{0.5 < |x| < 1\}}}

\code{parzen}:
\deqn{k(x) = (1 - 6x^2 + 6|x|^3) * 1_{\{0 \leq |x| \leq 0.5\}} + 2(1 - |x|)^3 * 1_{\{0.5 < |x| \leq 1\}}}

\code{QS}:
\deqn{k(x) = \frac{25}{12 \pi ^2 x^2} \left(\frac{\sin(6\pi x / 5)}{6\pi x / 5} - \cos(6 \pi x / 5)\right)}

\code{TH}:
\deqn{k(x) = (1 + \cos(\pi  x)) / 2 * 1_{\{|x| < 1\}}}

\code{truncated}:
\deqn{k(x) = 1_{\{|x| < 1\}}}

\code{SFT}:
\deqn{k(x) = (1 - 4(|x| - 0.5)^2)^2 * 1_{\{|x| < 1\}}}

\code{Epanechnikov}:
\deqn{k(x) =  3 \frac{1 - x^2}{4} * 1_{\{|x| < 1\}}}

\code{quatratic}:
\deqn{k(x) = (1 - x^2)^2 * 1_{\{|x| < 1\}}}
}

\value{
long run variance \eqn{\sigma^2} (numeric) resp. \eqn{\Sigma} (numeric matrix)
}

\author{
Sheila \enc{Görz}{Goerz}
}

\references{
Andrews, D.W. "Heteroskedasticity and autocorrelation consistent covariance matrix estimation." Econometrica: Journal of the Econometric Society (1991): 817-858.

Dehling, H., et al. "Change-point detection under dependence based on two-sample U-statistics." Asymptotic laws and methods in stochastics. Springer, New York, NY, (2015). 195-220.

Dehling, H., Fried, R., and Wendler, M. "A robust method for shift detection in time series." Biometrika 107.3 (2020): 647-660.

Parzen, E. "On consistent estimates of the spectrum of a stationary time series." The Annals of Mathematical Statistics (1957): 329-348.

Shao, X. "The dependent wild bootstrap." Journal of the American Statistical Association 105.489 (2010): 218-235.
}

\seealso{
\code{\link{CUSUM}}, \code{\link{HodgesLehmann}}, \code{\link{wilcox_stat}}
}

\examples{
Z <- c(rnorm(20), rnorm(20, 2))

## kernel density estimation
lrv(Z)

## overlapping subsampling
lrv(Z, method = "subsampling", control = list(overlapping = FALSE, distr = TRUE, l = 5))

## dependent wild bootstrap
lrv(Z, method = "bootstrap", control = list(B = 2000, l = 4, kFun = "parzen"))
}
