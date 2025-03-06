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

wilcox_stat <- function(x, h = 1L, method = "kernel", control = list())
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
  method <- match.arg(method, c("subsampling", "kernel", "bootstrap", "none"))
  
  n <- length(x)
  
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
      res <- CUSUM_cpp(x)
      # res <- .Call("CUSUM", as.numeric(x))
    } else
    {
      res <- .Call("wilcox", as.numeric(x))
    }
  } else if(is.function(h))
  {
    hVec <- Vectorize(h)
    res <- c(sum(hVec(x[1], x[-1])), rep(0, n-2))
    
    sapply(2:(n-1), function(i)
    {
      res[i] <<- res[i-1] - sum(hVec(x[1:(i-1)], x[i])) + sum(hVec(x[i], x[(i+1):n]))
    })
    res <- abs(res) / sqrt(n^3)
  } else
  {
    stop("Invalid argument for h!")
  }
  
  k <- which.max(res)
  
  if((method == "subsampling" & (is.null(control$l_n) || is.na(control$l_n))) | 
     (method == "kernel" & (is.null(control$b_n) || is.na(control$b_n))) | 
     (method == "bootstrap" & (is.null(control$l_n) || is.na(control$l_n))))
  {
    n <- length(x)
    x.adj <- x
    x.adj[(k+1):n] <- x.adj[(k+1):n] - mean(x[(k+1):n]) + mean(x[1:k])
    rho <- cor(x.adj[-n], x.adj[-1], method = "spearman")
    
    #####
    p1 <- ifelse(is.numeric(h) && h == 1, 0.25, 0.4)
    p2 <- ifelse(is.numeric(h) && h == 1, 0.8, 1/3)
    #####
    
    param <- max(ceiling(n^(p1) * ((2 * rho) / (1 - rho^2))^(p2)), 1)
    param <- min(param, n-1)
    if(is.na(param)) param <- 1
    control$b_n <- param
    control$l_n <- control$b_n
  }
  
  sigma <- sqrt(lrv(x, method = method, control = control))
  Tn <- max(res) / sigma
    
  attr(Tn, "cp-location") <- k
  attr(Tn, "data") <- ts(x)
  attr(Tn, "lrv-method") <- method
  attr(Tn, "sigma") <- sigma
  if(method == "kernel") attr(Tn, "param") <- control$b_n else
    attr(Tn, "param") <- control$l_n
  attr(Tn, "teststat") <- ts(res / sigma)
  class(Tn) <- "cpStat"

  return(Tn)
}





