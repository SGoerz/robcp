shift_cusum <- function(x, version = "empVar", alpha = 0.5, method = "kernel",
                        control = list(), fpc = TRUE, tol = 1e-8)
{
  Dataname <- deparse(substitute(x))
  
  
  if(version == "empVar")
  {
    stat <- CUSUM(x^2, method = method, control = control)
  } else 
  {
    n <- length(n)
    
    if(version == "MD")
    {
      require(cumstats)
      y <- cummedian(x)
      res <- .Call("MD", as.numeric(x), as.numeric(y), as.numeric(n))
      stat <- res / (1:n) - res[n] / n
    } else if(version == "GMD")
    {
      res <- .Call("GMD", as.numeric(x), as.numeric(n))
      stat <- res / ((1:n)^2) - res / (n^2)
    } else if(version == "MAD")
    {
      y <- sapply(1:n, function(k) mad(x[1:k]))
      stat <- y - y[n]
    } else if(version == "QAlpha")
    {
      sapply(2:(n-1), function(k)
      {
        a <- floor(k * (k - 1) / 2 * (1 - alpha)) + 1
        kthPair(x[1:k], -x[1:k], a)
      })
    } else
    {
      stop("Supplied version not supported.")
    }
    
    stat <- max(1:n * abs(stat)) / sqrt(n)
  } 
  
  location <- attr(stat, "cp-location")
  names(stat) <- "S"
  
  if(fpc) stat <- stat + 1.46035 / sqrt(2 * pi) / sqrt(length(x))
  
  erg <- list(alternative = "two-sided", method = "Huberized CUSUM test",
              data.name = Dataname, statistic = stat,
              p.value = 1 - pKSdist(stat, tol), 
              cp.location = location)
  
  class(erg) <- "htest"
  return(erg)
}


