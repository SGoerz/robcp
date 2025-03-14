##'pKSdist: asymptotic cumulative distribution of the maximal absolute value of
##'         of a Brownian bridge
##' 
##'input: tn (test statistic; numeric vector)
##'       tol (tolerance; numeric > 0)
##'       
##'output: vector of P(t_n(X) <= tn)

pKSdist <- function(tn, tol = 1e-8)
{
  if(!is.numeric(tol) || tol <= 0) stop("tol has to be a positive numeric!")
  erg <- .Call("pKSdist", as.numeric(tn), as.numeric(tol))
  return(erg)
}