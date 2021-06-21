## Hodges-Lehmann

##'@name HodgesLehmann
##'@title Hodges Lehmann Test Statistic
##'@description Computes the test statistic for the Hodges-Lehmann change point test.
##'@param x time series (numeric or ts vector).
##'@param b bandwidth for u_hat() (numeric > 0).
##'@param method method for estimating the long run variance.
##'@param control a list of control parameters.
##'@return Test statistic (numeric value) with the attribute cp-location 
##'        indicating at which index a change point is most likely. Is an S3 
##'        object of the class cpStat        
HodgesLehmann <- function(x, method = "subsampling", control = list())
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
  ## end argument check
  
  if(is.null(control$b2))
  {
    b <- control$b_n
  } else
  {
    b <- control$b2
    control$b2 <- NULL
  }
  
  n <- length(x)
  Mn <- sapply(1:(n-1), function(k)
  {
    medDiff <- medianDiff(x[(k+1):n], x[1:k])
    u_hat(x - c(rep(0, k), rep(medDiff, n - k)), b) *
      k / n * (1 - k / n) * abs(medDiff) 
  })
  Tn <- sqrt(n) * max(Mn) / sqrt(lrv(x, method = method, control = control))
  
  attr(Tn, "cp-location") <- which.max(Mn)
  class(Tn) <- "cpStat"
  
  return(Tn)
}


## default values?
u_hat <- function(x, b, kFun = "bartlett")
{
  if(b <= 0) 
    stop("b must be numeric, greater than 0 and smaller than the length of the time series!")
  
  n <- length(x)
  kFun <- pmatch(kFun, c("bartlett", "FT", "parzen", "QS", "TH", "truncated"))
  res <- .Call("u_hat", as.numeric(x), as.numeric(b), as.numeric(kFun))
  return(res)
}