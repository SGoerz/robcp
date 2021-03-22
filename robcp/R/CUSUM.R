##'teststat: computes the test statistic for the CUSUM test
##'
##'input: y (time series; vector, matrix or ts object)
##'       fun (psi function; character string)
##'       b_n (bandwidth for the long run variance estimation; numeric > 0)


CUSUM <- function(x, fun = "HLm", b_n, inverse = "Cholesky", ...)
{
  ## argument check
  if(is(x, "ts"))
  {
    class(x) <- "numeric"
  }
  if(!(is(x, "matrix") || is(x, "numeric") || is(x, "integer")))
  {
    stop("x must be a numeric or integer vector or matrix!")
  }
  ## end argument check
  
  #x <- psi(y, fun, k, constant)
  
  if(is(x, "matrix"))
  {
    sigma <- lrv(x, b_n, ...)
    
    m <- ncol(x)
    n <- nrow(x)
    
    #browser()
    
    if(inverse == "Cholesky")
    {
      # Modified Cholesky Factorization
      mchol <- modifChol(sigma, ...)
      swaps <- attr(mchol, "swaps")
      # invert sigma
      sigma.inv <- chol2inv(mchol) 
    } else if(inverse == "generalized")
    {
      if(!requireNamespace("MASS", quietly = TRUE)) 
      {
        stop("Package \"MASS\" needed for this function to work. Please install it.",
             call. = FALSE)
      }
      sigma.inv <- ginv(sigma, ...)
      swaps <- 0:(m - 1)
    }
    
    erg <- .Call("CUSUM_ma", as.numeric(x), as.numeric(sigma.inv), 
                 as.numeric(swaps), as.numeric(n), as.numeric(m))
    
    return(erg)
  }
  
  sigma.inv <- sqrt(lrv(x, b_n))
  erg <- .Call("CUSUM", as.numeric(x), as.numeric(sigma.inv))
  return(erg)
}