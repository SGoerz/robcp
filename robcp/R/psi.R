# HLm: marginal huberized location
# HLg: global huberized location
# SLm: marginal sign location
# SLg: global sign location
# HCm: marginal huberized covariance
# HCg: global huberized covariance
# SCm: marginal sign covariance
# SCg: global sign covariance

## psi: standartizes and transforms a time series Y according to a psi 
##      function 
##              
## input: Y (time series; numeric vector, matrix or ts object)
##        fun (psi function; charater, one of the above)
##        k (parameter for psi function; numeric)
##        constant (scale factor of the MAD; numeric)
##      
## output: transformed time series (numeric vector or matrix)
psi <- function(y, fun = c("HLm", "HLg", "SLm", "SLg", "HCm", "HCg", "SCm", "SCg"), 
                k, constant = 1.4826)
{
  ## argument check
  ### preserve (dim)names????
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
  
  fun <- match.arg(fun)
  fun <- which(c("HLm", "HLg", "SLm", "SLg", "HCm", "HCg", "SCm", "SCg") == fun)
  if(length(fun) == 0) stop("Wrong argument to \'fun\'")
  
  if(missing(k) || !is(k, "numeric"))
  {
    if(fun == 2 || fun == 6)
    {
      k <- sqrt(qchisq(0.8, df = ifelse(is(y, "matrix"), ncol(y), 1)))
    }
    else k <- 1.5
  }
  ## end argument check
  
  if(is(y, "matrix"))
  {
    m <- ncol(y)
    n <- nrow(y)
    med <- apply(y, 2, median)
    MAD <- apply(y, 2, function(x) mad(x, constant = constant))
  } else
  {
    if(fun %in% 7:8) 
    {
      stop("SCm and SCg are not available for a one-dimensional time series!")
    }
    n <- length(y)
    m <- 1
    med <- median(y)
    MAD <- mad(y, constant = constant)
  }
  
  erg <- .Call("psi", as.numeric(y), as.numeric(fun), as.numeric(n), 
               as.numeric(m), as.numeric(k), #as.numeric(constant), 
               # median and mad are being computed in R 
               as.numeric(med), as.numeric(MAD))

  if(is(y, "matrix")) erg <- matrix(erg, nrow = n)
  if(timeseries) 
  {
    erg <- ts(erg, start = tsp[1], end = tsp[2], frequency = tsp[3])
  }
  
  return(erg)
}