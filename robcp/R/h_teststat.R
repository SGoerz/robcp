##'h_teststat: computes the test statistic for the Huberized change point test
##'
##'input: y (time series; vector, matrix or ts object)
##'       fun (psi function; character string)
##'       b_n (bandwidth for the long run variance estimation; numeric > 0)
##'       k (parameter for the psi function; numeric > 0)
##'       constant (scaling factor of the mad; numerc > 0)


h_teststat <- function(y, fun = "HLm", b_n, k, constant = 1.4826)
{
  ## argument check
  if(is(y, "ts"))
  {
    class(y) <- "numeric"
  }
  if(!(is(y, "matrix") || is(y, "numeric") || is(y, "integer")))
  {
    stop("x must be a numeric or integer vector or matrix!")
  }
  ## end argument check
  
  x <- psi(y, fun, k, constant)
  
  if(is(x, "matrix"))
  {
    sigma <- lrv(x, b_n)
    
    # Modified Cholesky Factorization
    mchol <- modifChol(sigma)
    swaps <- attr(mchol, "swaps")
    # invert sigma
    sigma <- chol2inv(mchol)  
    
    m <- ncol(x)
    n <- nrow(x)
    
    erg <- .Call("h_teststat_ma", as.numeric(x), as.numeric(sigma), 
                 as.numeric(swaps), as.numeric(n), as.numeric(m))
    
    return(erg)
  }
  
  sigma <- sqrt(lrv(x, b_n))
  erg <- .Call("h_teststat", as.numeric(x), as.numeric(sigma))
  return(erg)
}