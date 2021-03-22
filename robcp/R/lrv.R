##'lrv: estimates the long run variance resp. covariance matrix of the supplied
##'     time series
##'
##'input: x (time series; numeric vector, matrix or ts object)
##'       b_n (bandwidth; numeric; default: length of time series ^ (1/3)
##'       
##'output: long run variance (numeric value) or long run covariance matrix 
##'        (numeric matrix with dim. m x m, when m is the number of columns)

lrv <- function(x, b_n, gamma0 = TRUE)
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
  if(!missing(b_n) && (!is(b_n, "numeric") || b_n <= 0))
  {
    stop("b_n must be numeric and greater than 0!")
  }
  ## end argument check
  
  if(is(x, "matrix"))
  {
    m <- ncol(x)
    n <- nrow(x)

    if(missing(b_n)) b_n <- n^(1/3)
    if(b_n > n)
      stop("The bandwidth b_n cannot be larger than the length of the time series!")
    
    x_cen <- apply(x, 2, function(x) x - mean(x))
    
    erg <- .Call("lrv_matrix", as.numeric(x_cen), 
                 as.numeric(n), as.numeric(m), as.numeric(b_n), 
                 PACKAGE = "robcp")
    
    erg <- matrix(erg, ncol = m)
  } else
  {
    n <- length(x)
    if(missing(b_n)) b_n <- n^(1/3)
    if(b_n > n) 
      stop("The bandwidth b_n cannot be larger than the length of the time series!")
    
    x_cen <- x - mean(x)
    erg <- .Call("lrv", as.numeric(x_cen), as.numeric(b_n),
                 PACKAGE = "robcp")
    if(erg < 0 & gamma0)
    {
      warning("Estimated long run variance was < 0; only the estimated autocovariance to lag 0 is returned!")
      erg <- (n - 1) / n * var(x)
    }
  }
  
  return(erg)
}