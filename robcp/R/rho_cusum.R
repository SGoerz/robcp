## Test for fluctuation changes based on CUSUM and Spearman's rho

##'@name scale_stat
##'@title Test statistic to detect changes in Spearman's rho
##'@description Computes the test statistic for a CUSUM-based test on fluctuation changes.
##'@param x time series (numeric or ts vector).
##'@param control a list of control parameters.
##'@return Test statistic (numeric value) with the attribute cp-location 
##'        indicating at which index a change point is most likely. Is an S3 
##'        object of the class cpStat   

rho_stat <- function(x, control = list())
{
  ## argument check
  if(is(x, "ts"))
  {
    class(x) <- "numeric"
    tsp <- attr(x, "tsp")
  }
  
  ######## also for 1-dim data? ########
  # if(!(is(x, "matrix") || is(x, "numeric") || is(x, "integer")))
  # {
  #   stop("x must be a numeric or integer vector or matrix!")
  # }
  # if(length(x) < 2)
  # {
  #   stop("x must consist of at least 2 observations!")
  # }
  # 
  # 
  if(is.null(control$kFun) || is.na(control$kFun))
  {
    control$kFun <- "bartlett"
  }
  n <- nrow(x)
  if(is.null(control$b_n) || is.na(control$b_n))
  {
    control$b_n <- log(n)
  }
  
  ## end argument check
  
  d <- ncol(x)
  rks <- 1 - apply(x, 2, rank) / n
  prd <- apply(rks, 1, prod)
  res <- (2^d / (1:n) * cumsum(prd) - 1) * (d + 1) / (2^d - d - 1)
      
  control$version <- "rho"

  stat <- res - res[n]
  stat <- 1:n * abs(stat)
  k <- which.max(stat) + 1
  res <- max(stat) / sqrt(n * lrv(x, method = "kernel", control = control))
   
  attr(res, "cp-location") <- as.integer(k)
  class(res) <- "cpStat"
  
  return(res)
}


##'@name rho_cusum
##'@title A fluctuation test for constant Spearman’s rho with nuisance-free limit distribution
##'@description Performs a CUSUM-based test on changes in Spearman's rho.
##'@param x time series (numeric or ts vector).
##'@param control a list of control parameters.
##'@param constant scale factor for the MAD.
##'@param fpc finite population correction (boolean).
##'@param tol tolerance of the distribution function (numeric), which is used do compute p-values.
##'@return A list fo the class "htest" containing

rho_cusum <- function(x, control = list(), fpc = TRUE, tol = 1e-8)
{
  Dataname <- deparse(substitute(x))
  
  stat <- rho_stat(x = x, control = control)
  location <- attr(stat, "cp-location")
  names(stat) <- "S"
  
  if(fpc) stat <- stat + 1.46035 / sqrt(2 * pi) / sqrt(length(x))
  
  erg <- list(alternative = "two-sided", method = "CUSUM test for changes in Spearman's rho",
              data.name = Dataname, statistic = stat,
              p.value = 1 - pKSdist(stat, tol), 
              cp.location = location)
  
  class(erg) <- "htest"
  return(erg)
}
