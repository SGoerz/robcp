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
HodgesLehmann <- function(x, b_u = NA, method = "subsampling", control = list())
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
  n <- length(x)
  
  ## first iteration (k == 1)
  medDiff <- medianDiff(x[2:n], x[1])
  x.adj <- x - c(0, rep(medDiff, n - 1))
  
  ## determine b_u adaptively
  if(is.na(b_u))
  {
    # diffs <- unlist(sapply(1:(n-1), function(i)
    # {
    #   x.adj[(i+1):n] - x.adj[i]
    # }))
    
    diffs <- unlist(sapply(1:(n-1), function(i)
    {
      temp <- rep(x.adj[(i+1):n], each = i)
      x.adj[1:i] - temp
    }))
    
    ## sometimes diffs is too large
    b_u <- tryCatch(bw.SJ(diffs), error = function(e) bw.nrd0(diffs))
  }
  
  ## first Mn
  Mn <- u_hat(x.adj, b_u, "QS") * (n-1) / n^2 * abs(medDiff)
  
  ## next Mn's
  Mn <- c(Mn, sapply(2:(n-1), function(k)
  {
    medDiff <- medianDiff(x[(k+1):n], x[1:k])
    x.adj <- x - c(rep(0, k), rep(medDiff, n - k))
    
    u_hat(x.adj, b_u, "QS") *
      k / n * (1 - k / n) * abs(medDiff) 
  }))
  
  k <- which.max(Mn)
  
  if(method == "subsampling" & (is.null(control$l) || is.na(control$l)))
  {
    x.adj <- x
    x.adj[(k+1):n] <- x.adj[(k+1):n] - mean(x[(k+1):n]) + mean(x[1:k])
    rho <- cor(x.adj[-n], x.adj[-1], method = "spearman")
    control$l <- max(ceiling(n^(1/3) * ((2 * rho) / (1 - rho^2))^(2/3)), 1)
  }
  
  Tn <- sqrt(n) * max(Mn) / sqrt(lrv(x, method = method, control = control))
  
  attr(Tn, "cp-location") <- k
  class(Tn) <- "cpStat"
  
  return(Tn)
}



u_hat <- function(x, b_u, kFun = "QS")
{
  if(b_u <= 0) 
    stop("b must be numeric, greater than 0 and smaller than the length of the time series!")
  
  n <- length(x)
  kFun <- pmatch(kFun, c("bartlett", "FT", "parzen", "QS", "TH", "truncated",
                         "Gaussian"))
  res <- .Call("u_hat", as.numeric(x), as.numeric(b_u), as.numeric(kFun))
  return(res)
}