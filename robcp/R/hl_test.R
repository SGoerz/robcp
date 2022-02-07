## hl_test

##'@name hl_test
##'@title Hodges-Lehmann Test for Change Points 
##'@description Performs the Hodges-Lehmann change point test.
##'@param x time series (numeric or ts vector).
##'@param method method for estimating the long run variance.
##'@param control a list of control parameters.
##'@param tol tolerance of the distribution function (numeric), which is used do compute p-values.
##'@return A list fo the class "htest" containing
hl_test <- function(x, b_u = "nrd0", method = "subsampling", control = list(), 
                    tol = 1e-8, plot = FALSE)
{
  Dataname <- deparse(substitute(x))
  stat <- HodgesLehmann(x, b_u = b_u, method = method, control = control)
  location <- attr(stat, "cp-location")
  names(stat) <- "S"
  
  erg <- list(alternative = "two-sided", method = "Hodges-Lehmann change point test",
              data.name = Dataname, statistic = stat,
              p.value = 1 - pKSdist(stat, tol), 
              cp.location = location, 
              lrv = list(method = attr(stat, "lrv-method"), 
                         param = attr(stat, "param"), 
                         value = attr(stat, "sigma")))
  
  if(plot) plot(stat)
  
  class(erg) <- "htest"
  return(erg)
}