##'lrv: estimates the long run variance resp. covariance matrix of the supplied
##'     time series
##'
##'input: x (time series; numeric vector, matrix or ts object)
##'       b_n (bandwidth; numeric; default: length of time series ^ (1/3)
##'       
##'output: long run variance (numeric value) or long run covariance matrix 
##'        (numeric matrix with dim. m x m, when m is the number of columns)

lrv <- function(x, b_n)
{
  ## argument check
  if(is(x, "ts"))
  {
    timeseries <- TRUE
    tsp <- attr(x, "tsp")
    class(x) <- "numeric"
  } else
  {
    timeseries <- FALSE
    tsp <- NULL
  }
  if(!(is(x, "matrix") || is(x, "numeric") || is(x, "integer")))
  {
    stop("x must be a numeric or integer vector or matrix!")
  }
  if(!missing(b_n) && (!is(bn, "numeric") || b_n <= 0))
  {
    stop("b_n must be numeric and greater than 0!")
  }
  ## end argument check
  
  if(is(x, "matrix"))
  {
    m <- ncol(x)
    n <- nrow(x)
    
    if(missing(b_n)) b_n <- n^(1/3)
    
    x_cen <- apply(x, 2, function(x) x - mean(x))
    
    erg <- .Call("lrv_matrix", as.numeric(x_cen), 
                 as.numeric(n), as.numeric(m), as.numeric(b_n), 
                 PACKAGE = "robcp")
    
    erg <- matrix(erg, ncol = m)
  } else
  {
    if(missing(b_n)) b_n <- length(x)^(1/3)
    x_cen <- x - mean(x)
    erg <- .Call("lrv", as.numeric(x_cen), as.numeric(b_n),
                 PACKAGE = "robcp")
  }
  
  if(timeseries)
  {
    erg <- ts(erg, start = tsp[1], end = tsp[2], frequency = tsp[3])
  }
  
  return(erg)
}