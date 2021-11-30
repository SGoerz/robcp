## Tests for Scale Changes Based on Pairwise Differences (Gerstenberger et. al.)

##'@name scale_stat
##'@title Test statistic to detect Scale Changes 
##'@description Computes the test statistic for CUSUM-based tests on scale changes.
##'@param x time series (numeric or ts vector).
##'@param version variance estimation method. One of "empVar", "MD", "GMD", "MAD", "QBeta".
##'@param method method for estimating the long run variance.
##'@param control a list of control parameters.
##'@param constant scale factor for the MAD.
##'@param beta quantile of the distribution function of all absolute pairwise differences used in \code{version = "QBeta"}.
##'@return Test statistic (numeric value) with the attribute cp-location 
##'        indicating at which index a change point is most likely. Is an S3 
##'        object of the class cpStat   

scale_stat <- function(x, version = "empVar", method = "kernel", 
                       control = list(), constant = 1.4826, beta = 0.5)
{
  ## argument check
  if(is(x, "ts"))
  {
    class(x) <- "numeric"
    tsp <- attr(x, "tsp")
  }
  if(!(is(x, "matrix") || is(x, "numeric") || is(x, "integer")))
  {
    stop("x must be a numeric or integer vector or matrix!")
  }
  if(length(x) < 2)
  {
    stop("x must consist of at least 2 observations!")
  }
  ## end argument check
  
  
  if(version == "empVar")
  {
    return(CUSUM(x^2, method = method, control = control))
  } else 
  {
    n <- length(x)
    
    if(version == "MD")
    {
      require(cumstats)
      y <- cummedian(x)
      res <- .Call("MD", as.numeric(x), as.numeric(y), as.numeric(n))
      stat <- res / (2:n) - res[n-1] / n
    } else if(version == "GMD")
    {
      res <- .Call("GMD", as.numeric(x), as.numeric(n))
      stat <- res / ((2:n)^2) - res[n-1] / (n^2)
    } else if(version == "MAD")
    {
      y <- sapply(2:n, function(k) mad(x[1:k]), constant = constant)
      stat <- y - y[n-1]
    } else if(version == "QAlpha")
    {
      stat <- sapply(1:(n-1), function(k)
      {
        a <- floor(k * (n - k) * (1 - beta))
        kthPair(x[1:k], -x[1:k], a)
      })
    } else
    {
      stop("Supplied version not supported.")
    }
    
    stat <- 2:n * abs(stat)
    print(stat)
    k <- which.max(stat) + 1
    res <- max(2:n * abs(stat)) / sqrt(n) / lrv(x, method = method, control = control)
    
    attr(res, "cp-location") <- as.integer(k)
    class(res) <- "cpStat"
    
    return(res)
  } 
}

## Scale change test

##'@name scale_cusum
##'@title Tests for Scale Changes Based on Pairwise Differences
##'@description Performs the CUSUM-based test on changes in the scale.
##'@param x time series (numeric or ts vector).
##'@param version variance estimation method. One of "empVar", "MD", "GMD", "MAD", "QBeta".
##'@param method method for estimating the long run variance.
##'@param control a list of control parameters.
##'@param constant scale factor for the MAD.
##'@param beta quantile of the distribution function of all absolute pairwise differences used in \code{version = "QBeta"}.
##'@param fpc finite population correction (boolean).
##'@param tol tolerance of the distribution function (numeric), which is used do compute p-values.
##'@return A list fo the class "htest" containing

scale_cusum <- function(x, version = "empVar", method = "kernel",
                        control = list(), constant = 1.4826, beta = 0.5, 
                        fpc = TRUE, tol = 1e-8)
{
  Dataname <- deparse(substitute(x))
  
  stat <- scale_stat(x = x, version = version, method = method,
                     control = control, constant = constant, beta = beta)
  location <- attr(stat, "cp-location")
  names(stat) <- "S"
  
  if(fpc) stat <- stat + 1.46035 / sqrt(2 * pi) / sqrt(length(x))
  
  erg <- list(alternative = "two-sided", method = "CUSUM test for scale changes",
              data.name = Dataname, statistic = stat,
              p.value = 1 - pKSdist(stat, tol), 
              cp.location = location)
  
  class(erg) <- "htest"
  return(erg)
}

f <- function(n, alpha) 
{
  sapply(1:n, function(k) floor(k * (n - k) * (1 - alpha)))
}

f(5, 0.5)
