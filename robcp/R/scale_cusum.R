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

scale_stat <- function(x, version = "empVar", control = list(), 
                       constant = 1.4826, beta = 0.5)
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
  if(is.null(control$kFun) || is.na(control$kFun))
  {
    control$kFun <- "SFT"
  }
  n <- length(x)
  if(is.null(control$b_n) || is.na(control$b_n))
  {
    rho <- abs(acf(x, plot = FALSE)[[1]][, , 1])
    kappa <- max(5, sqrt(log10(n)))
    i <- 1
    cond <- 2 * sqrt(log10(n) / n)
    repeat
    {
      if(max(rho[i:(kappa+i)]) <= cond) break
      i <- i + 1
    }
    control$b_n <- i - 1
  }
  
  ## end argument check
  
  if(version == "empVar")
  {
    control$mean <- mean(x)
    control$var <- var(x)
    
    temp <- .Call("CUSUM", as.numeric(x^2))
    res <- temp[1] / sqrt(lrv(x, method = "scale_kernel", control = control))
    k <- temp[2]
  } else 
  {
    if(version == "MD")
    {
      require(cumstats)
      y <- cummedian(x)
      res <- .Call("MD", as.numeric(x), as.numeric(y), as.numeric(n)) / (1:(n-1))

      control$mean <- median(x)
    } else if(version == "GMD")
    {
      res <- .Call("GMD", as.numeric(x), as.numeric(n)) / ((1:(n-1)) * (2:n)) * 2
    } else if(version == "MAD")
    {
      res <- sapply(2:n, function(k) mad(x[1:k], constant = constant))
      control$mean <- median(x)
    } else if(version == "QBeta")
    {
      res <- QBeta(x, beta)
      # sorted <- x[1]
      # res <- sapply(2:n, function(k)
      # {
      #   sorted <<- sort(c(sorted, x[k]), decreasing = TRUE)
      #   a <- ceiling(k * (k - 1) / 2 * (1 - beta))
      #   # a <- floor(k * (k - 1) / 2 * (1 - beta)) + 1
      #   kthPair(sorted[1:(k-1)], -sorted[2:k], a)
      # })
      control$mean <- beta
    } else
    {
      stop("version not supported.")
    }
    
    control$var <- res[n-1]
    control$version <- version
    
    stat <- res - res[n-1]
    stat <- 2:n * abs(stat)
    k <- which.max(stat) + 1
    res <- max(stat) / sqrt(n * lrv(x, method = "scale_kernel", control = control))
  } 
  attr(res, "cp-location") <- as.integer(k)
  class(res) <- "cpStat"
  
  return(res)
}

f <- function(k, a)
{
  ceiling(k * (k - 1) / 2 * (1 - a))
}

f(2, 0.1)


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

scale_cusum <- function(x, version = "empVar", control = list(), 
                        constant = 1.4826, beta = 0.5, fpc = TRUE, tol = 1e-8)
{
  Dataname <- deparse(substitute(x))
  
  stat <- scale_stat(x = x, version = version, control = control,
                     constant = constant, beta = beta)
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
