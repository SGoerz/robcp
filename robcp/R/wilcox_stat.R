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
  if(is.numeric(h))
  {
    if(!(h %in% 1:2)) 
    {
      stop("Wrong test version h!")
      #stop("h must be either 1 or 2!")
    }
    if(is.null(control$distr))
    {
      control$distr <- h == 1L
    } else
    {
      if(is.logical(control$distr) & 2 - as.numeric(control$distr) != h)
      {
        warning("Argument for 'distr' not suitable for test version.")
      }
    }
    if(h == 2L)
    {
      res <- .Call("CUSUM", as.numeric(x))
    } else
    {
      res <- .Call("wilcox", as.numeric(x), as.numeric(h))
    }
  } else if(is.function(h))
  {
    n <- length(x)
    hVec <- Vectorize(h)
    res <- sum(hVec(x[1], x[-1]))
    max <- abs(res)
    loc <- 1
    
    sapply(2:(n-1), function(k)
    {
      res <<- res - sum(hVec(x[1:(k-1)], x[k])) + sum(hVec(x[k], x[(k+1):n]))
      
      if(abs(res) > max)
      {
        max <<- abs(res)
        loc <<- k
      }
    })
    res[1] <- max
    res[2] <- loc
  } else
  {
    stop("Invalid argument for h!")
  }
  
  Tn <- res[1] / sqrt(lrv(x, method = method, control = control))
    
  attr(Tn, "cp-location") <- res[2]
  class(Tn) <- "cpStat"
  
  return(Tn)
}
