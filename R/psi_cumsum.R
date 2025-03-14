## psi_cumsum: cumulative sums of psi-transformed (multivariate) time series
##
## input: y (time series; numeric vector, matrix or ts object)
##        fun (psi function; character)
##        k (parameter for psi function; numeric)
##        constant (scale factor of the MAD; numeric)
##
## output: cumulative sum vector, matrix or ts object
psi_cumsum <- function(y, fun = "HLm", k, constant = 1.4826)
{
  ## argument check
  if(is(y, "ts"))
  {
    timeseries <- TRUE
    tsp <- attr(y, "tsp")
    class(y) <- "numeric"
  } else
  {
    timeseries <- FALSE
    tsp <- NULL
  }
  if(!(is(y, "matrix") | is(y, "numeric") | is(y, "integer")))
  {
    stop("x must be numeric and either be a vector, a matrix or a ts object!")
  }
  ## end argument check
  
  x <- psi(y, fun, k, constant)
  
  if(is(x, "matrix"))
  {
    m <- ncol(x)
    n <- nrow(x)
    
    # erg <- cumsum_ma_cpp(x)
    erg <- .Call("c_cumsum_ma", as.numeric(x), as.numeric(n), as.numeric(m))
    erg <- matrix(erg, ncol = m)
  }
  else
  {
    erg <- cumsum(x)
  }
  
  if(timeseries)
  {
    erg <- ts(erg, start = tsp[1], end = tsp[2], frequency = tsp[3])
  }
  
  return(erg)
}