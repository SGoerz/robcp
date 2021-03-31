##'huber_cusum: performs a CUSUM test on data transformed (Huberized) by psi() 
##'
##'input: x (time series; numeric vector, matrix or ts object)
##'       fun (psi function; character)
##'       tol (tolerance of the distribution function)
##'       b_n (bandwidth for the long run variance computation; numeric > 0)
##'       k (parameter for the psi function; numeric > 0)
##'       constant (scaling factor for the mad; numeric > 0)
##'       fpc (finite population correction; boolean)
##'       tol (tolerance; numeric > 0)
##'       
##'output: htest object (list) containing:
##'        - alternative; method; data.name (character)
##'        - statistic (numeric)
##'        - p.value (numeric)

huber_cusum <- function(x, fun = "HLm", b_n, k, constant = 1.4826, fpc = TRUE, 
                        tol = 1e-8, ...)
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
  
  Dataname <- deparse(substitute(x))
  stat <- CUSUM(psi(x, fun = fun, k, constant), b_n, ...)
  location <- attr(stat, "cp-location")
  names(stat) <- "S"
  
  if(is(x, "matrix"))
  {
    h <- ncol(x)
    h <- switch(fun, HCm = h * (h + 1) / 2, HCg = h * (h + 1) / 2, 
                SCm = h * (h - 1) / 2, SCg = h * (h + 1) / 2 - 1, h)
    
    if(fpc) stat <- (sqrt(stat) + 1.46035 / sqrt(2 * pi) / sqrt(nrow(x)))^2
    
    erg <- list(alternative = "two-sided", method = "Huberized CUSUM test",
                data.name = Dataname, statistic = stat,
                p.value = 1 - pBessel(stat, h), 
                cp.location = location)
  }
  else
  {
    if(fpc) stat <- stat + 1.46035 / sqrt(2 * pi) / sqrt(length(x))
    
    erg <- list(alternative = "two-sided", method = "Huberized CUSUM test",
                data.name = Dataname, statistic = stat,
                p.value = 1 - pKSdist(stat, tol), 
                cp.location = location)
  }
  
  class(erg) <- "htest"
  return(erg)
}