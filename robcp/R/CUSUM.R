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
      sigma.inv <- MASS::ginv(sigma, ...)
      swaps <- 0:(m - 1)
    }
    
    temp <- .Call("CUSUM_ma", as.numeric(x), as.numeric(sigma.inv), 
                 as.numeric(swaps), as.numeric(n), as.numeric(m))
  } else
  {
    sigma.inv <- sqrt(lrv(x, method = method, control = control))
    temp <- .Call("CUSUM", as.numeric(x), as.numeric(sigma.inv))
  }
  
  erg <- temp[1]
  attr(erg, "cp-location") <- as.integer(temp[2])
  
  class(erg) <- "cpStat"
  
  return(erg)
}

##'print.cpStat: print method for change point statistics
##'              prints the value of the test statistic and add the most likely
##'              change point location
##'@name CUSUM
print.cpStat <- function(x, ...)
{
  loc <- attr(x, "cp-location")
  print(round(as.numeric(x), digits = getOption("digits")), ...)
  cat("location: ", loc)
  return(invisible(x))
}

