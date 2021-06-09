## Wilcoxon-Mann-Whitney 

##'@name wilcox_stat
##'@title Wilcoxon-Mann-Whitney Test Statistic for Change Points
##'@description Computes the test statistic for the Wilcoxon-Mann-Whitney change point test
##'@param x time series (numeric or ts vector)
##'@param h version of the test (integer, 1 or 2)
##'@param method method for estimating the long run variance
##'@param control a list of control parameters for the estimation of the long run variance
##'@return Test statistic (numeric value) with the attribute cp-location 
##'        indicating at which index a change point is most likely. Is an S3
##'        object of the class cpStat

wilcox_stat <- function(x, h = 1L, method = "subsampling", control = list())
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
  
  if(is.null(control$overlapping)) 
  {
    control$overlapping <- FALSE
  }
  if(is.null(control$distr))
  {
    control$distr <- h == 1L
  }
  ## end argument check
  
  res <- .Call("wilcox", as.numeric(x), as.numeric(h))
  Tn <- res[1] / sqrt(lrv(x, method = method, control = control))
    
  attr(Tn, "cp-location") <- res[2]
  class(Tn) <- "cpStat"
  
  return(Tn)
}
