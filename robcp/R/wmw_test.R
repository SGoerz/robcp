## wmw_test

##'@name wmw_test
##'@title Wilocxon-Mann-Whitney Test for Change Points 
##'@description Performs the Wilcoxon-Mann-Whitney change point test.
##'@param x time series (numeric or ts vector).
##'@param h version of the test (integer, 1 or 2)
##'@param method method for estimating the long run variance.
##'@param control a list of control parameters.
##'@param tol tolerance of the distribution function (numeric), which is used do compute p-values.
##'@return A list fo the class "htest" containing
wmw_test <- function(x, h = 1L, method = "subsampling", control = list(), tol = 1e-8)
{
  # ## argument check
  # if(is(x, "ts"))
  # {
  #   class(x) <- "numeric"
  # }
  # if(!(is(x, "matrix") || is(x, "numeric") || is(x, "integer")))
  # {
  #   stop("x must be a numeric or integer vector or matrix!")
  # }
  # ## end argument check
  
  Dataname <- deparse(substitute(x))
  stat <- wilcox_stat(x, h = h, method = method, control = control)
  location <- attr(stat, "cp-location")
  names(stat) <- "S"
  
  erg <- list(alternative = "two-sided", method = "Wilcoxon-Mann-Whitney change point test",
              data.name = Dataname, statistic = stat,
              p.value = 1 - pKSdist(stat, tol), 
              cp.location = location)
  
  class(erg) <- "htest"
  return(erg)
}