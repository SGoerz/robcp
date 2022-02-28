## Test for fluctuation changes based on CUSUM and Spearman's rho

##'@name cor_stat
##'@title Test statistic to detect changes in the correlation in a time series.
##'@description Computes the test statistic for a CUSUM-based test on fluctuation changes.
##'@param x time series (numeric or ts vector).
##'@param control a list of control parameters.
##'@return Test statistic (numeric value) with the attribute cp-location 
##'        indicating at which index a change point is most likely. Is an S3 
##'        object of the class cpStat   

cor_stat <- function(x, version = "rho", control = list())
{
  ## argument check
  if(is(x, "ts"))
  {
    class(x) <- "numeric"
    tsp <- attr(x, "tsp")
  }
  
  if(!is(x, "matrix") || ncol(x) > 2)
  {
    stop("x must be a numeric or integer matrix with at least 2 columns!")
  }
  ## end argument check
  
  n <- nrow(x)
  
  if(version == "rho")
  {
    if(is.null(control$b_n) || is.na(control$b_n))
    {
      control$b_n <- log(n)
    }
    if(is.null(control$kFun) || is.na(control$kFun))
    {
      control$kFun <- "bartlett"
    }
    
    d <- ncol(x)
    rks <- 1 - apply(x, 2, rank) / n
    prd <- apply(rks, 1, prod)
    res <- (2^d / (1:n) * cumsum(prd) - 1) * (d + 1) / (2^d - d - 1)
  } else if(version == "tau")
  {
    if(is.null(control$b_n) || is.na(control$b_n))
    {
      control$b_n <- floor(2 * n^(1/3))
    }
    if(is.null(control$kFun) || is.na(control$kFun))
    {
      control$kFun <- "quadratic"
    }
    
    if(is.null(dim(x)) || ncol(x) != 2) stop("For the test on Kendall's tau, x must be bivariate.")
    res <- .Call("tau", as.numeric(x[, 1]), as.numeric(x[, 2]), as.numeric(n))
    control$scale <- res[n-1]
  } else
  {
    stop("version not supported.")
  }
  
  control$version <- version

  stat <- res - tail(res, 1)
  stat <- (n - length(res) + 1):n * abs(stat)
  k <- which.max(stat) + 1
  sigma <- sqrt(n * lrv(x, method = "kernel", control = control))
  res <- max(stat) / sigma
   
  attr(res, "cp-location") <- as.integer(k)
  attr(res, "data") <- ts(x)
  attr(res, "lrv-estimation") <- "kernel"
  attr(res, "sigma") <- sigma
  attr(res, "b_n") <- control$b_n 
  attr(res, "kFun") <- control$kFun
  attr(res, "teststat") <- ts(stat / sigma)
  class(res) <- "cpStat"
  
  return(res)
}


##'@name cor_cusum
##'@title Test to detect changes in the correlation in a time series.
##'@description Performs a CUSUM-based test on changes in Spearman's rho or Kendall's tau
##'@param x time series (numeric or ts vector).
##'@param control a list of control parameters.
##'@param fpc finite population correction (boolean).
##'@param tol tolerance of the distribution function (numeric), which is used do compute p-values.
##'@return A list fo the class "htest"
cor_cusum <- function(x, version = "rho", control = list(), fpc = TRUE, 
                      tol = 1e-8, plot = FALSE)
{
  Dataname <- deparse(substitute(x))
  
  stat <- cor_stat(x = x, version = version, control = control)
  location <- attr(stat, "cp-location")
  names(stat) <- "S"
  
  if(fpc) stat <- stat + 1.46035 / sqrt(2 * pi) / sqrt(length(x))
  
  erg <- list(alternative = "two-sided", method = "CUSUM test for changes in the correlation",
              data.name = Dataname, statistic = stat,
              p.value = 1 - pKSdist(stat, tol), 
              cp.location = location, 
              lrv = list(method = attr(stat, "lrv-method"), 
                         param = attr(stat, "param"), 
                         value = attr(stat, "sigma")))
  
  if(plot) plot(stat)
  
  class(erg) <- "htest"
  return(erg)
}
