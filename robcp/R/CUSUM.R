##'CUSUM: computes the test statistic for the CUSUM test
##'
##'input: y (time series; vector, matrix or ts object)
##'       fun (psi function; character string)
##'       b_n (bandwidth for the long run variance estimation; numeric > 0)
##'       inverse (character string specifying the method of inversion)
##'       ... (further arguments for the inverse-computing functions)
##'       
##'output: test statistic (numeric) with the attribute "cp-location" indicating 
##'        at which index a change point is most likely
##'        -> class "cpStat"

CUSUM <- function(x, method = "kernel", control = list(), inverse = "Cholesky", ...)
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
  method <- match.arg(method, c("subsampling", "kernel", "bootstrap", "none"))
  ## end argument check
  
  
  if(is(x, "matrix"))
  {
    sigma <- lrv(x, method = method, control = control)
    
    m <- ncol(x)
    n <- nrow(x)
    
    if(inverse == "Cholesky")
    {
      # Modified Cholesky Factorization
      mchol <- modifChol(sigma, ...)
      # invert sigma
      sigma.inv <- chol2inv(mchol) 
    } else if(inverse == "generalized")
    {
      if(!requireNamespace("MASS", quietly = TRUE)) 
      {
        stop("Package \"MASS\" needed for this function to work. Please install it.",
             call. = FALSE)
      }
      sigma.inv <- MASS::ginv(sigma, ...)
    } else if(inverse == "svd")
    {
      temp <- svd(sigma)
      sigma.inv <- temp$u %*% diag(1 / sqrt(temp$d)) %*% t(temp$v)
    }
    
    swaps <- 0:(m - 1)
    temp <- CUSUM_ma_cpp(x, sigma.inv)
    # temp <- .Call("CUSUM_ma", as.numeric(x), as.numeric(sigma.inv), 
    #               as.numeric(swaps), as.numeric(n), as.numeric(m))
    k <- which.max(temp)
  } else
  {
    temp <- CUSUM_cpp(x)
    # temp <- .Call("CUSUM", as.numeric(x))
    k <- which.max(temp)
    
    if((method == "subsampling" & (is.null(control$l_n) || is.na(control$l_n))) | 
       (method == "kernel" & (is.null(control$b_n) || is.na(control$b_n))) | 
       (method == "bootstrap" & (is.null(control$l_n) || is.na(control$l_n))))
    {
      n <- length(x)
      x.adj <- x
      x.adj[(k+1):n] <- x.adj[(k+1):n] - mean(x[(k+1):n]) + mean(x[1:k])
      rho <- abs(cor(x.adj[-n], x.adj[-1], method = "spearman"))

      ###
      p1 <- 0.45
      p2 <- 0.4
      ###
      
      param <- max(ceiling(n^(p1) * ((2 * rho) / (1 - rho^2))^(p2)), 1)
      param <- min(param, n-1)
      if(is.na(param)) param <- 1
      
      control$b_n <- param
      control$l_n <- param
    }
    
    if(method == "kernel" & (is.null(control$kFun) || is.na(control$kFun)))
    {
      control$kFun <- "TH"
    }
    
    sigma <- sqrt(lrv(x, method = method, control = control))
    temp <- temp / sigma
  }
  
  erg <- max(temp)
  attr(erg, "cp-location") <- k
  attr(erg, "data") <- ts(x)
  attr(erg, "lrv-method") <- method
  attr(erg, "sigma") <- sigma
  if(method == "kernel") attr(erg, "param") <- control$b_n else
    attr(erg, "param") <- control$l
  attr(erg, "teststat") <- ts(temp)
  if(is(x, "matrix")) attr(erg, "m") <- m
  
  class(erg) <- "cpStat"
  
  return(erg)
}


