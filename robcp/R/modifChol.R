##'modifChol: Computes the revised modified Cholesky factorization described in
##'           Schnabel and Eskow (1999)
##'           
##'input: x (symmetric numeric matrix)
##'       tau ()
##'       tau_bar ()
##'       mu ()
##'       
##'output: lower triangular matrix L of the form LL' = x + E with the attribute
##'        "swaps", a vector indicating which rows and columns were swapped 
##'        (more on help page)

modifChol <- function(x, tau = .Machine$double.eps^(1/3),
                      tau_bar = .Machine$double.eps^(2/3), mu = 0.1)
{
  ## argument check
  ### !! more for tau, tau_bar and mu !!
  if(!is(x, "matrix") | !is.numeric(x)) stop("x must be a numeric matrix!")
  
  n <- nrow(x)
  if(n != nrow(x)) stop("x must be a square matrix!")
  ## end argument check 
  
  erg <- matrix(.Call("cholesky", as.numeric(x), as.numeric(n),as.numeric(tau),
                      as.numeric(tau_bar), as.numeric(mu)), nrow = n)
  
  swaps <- erg[, n + 1]
  erg <- erg[, -(n + 1)]
  attr(erg, "swaps") <- swaps
  return(erg)
}