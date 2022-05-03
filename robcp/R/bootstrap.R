estim <- function(x, version, constant = 1.4826)#, alpha = 0.5)
{
  n <- length(x)
  
  if(version == "empVar")
  {
    return(max(.Call("CUSUM", as.numeric(x^2))))
  } else if(version == "MD")
  {
    y <- cumstats::cummedian(x)
    res <- .Call("MD", as.numeric(x), as.numeric(y), as.numeric(n)) / (1:(n-1))
  } else if(version == "GMD")
  {
    res <- .Call("GMD", as.numeric(x), as.numeric(n)) / ((1:(n-1)) * (2:n)) * 2
  # } else if(version == "MAD")
  # {
  #   res <- sapply(2:n, function(k) mad(x[1:k], constant = constant))
  # } else if(version == "Qalpha")
  # {
  #   res <- Qalpha(x, alpha)
  }
  return(max(abs(res - res[n-1]) * 2:n) / sqrt(n))
}


##'@description Performs a dependent block bootstrap in order to determine the p-value for the 
##'
##'@param x time series
##'@param l block length, 1 <= l
##'@param B number of bootstrap samples, numeric > 0
##'@param seed start for random number generator
##'       
##'@return 
dbb <- function(stat, data, version, control = list(), #alpha = 0.5, 
                constant = 1.4826, level = 0.05)
{
  n <- length(data)
  
  ## if l, B or seed are missing
  
  B <- control$B
  l <- control$l
  seed <- control$seed
  
  if(!is.numeric(B) || B < 1)
  {
    stop("B has to be a positive integer!")
  }
  if(!is.numeric(l) || l < 1)
  {
    stop("l has to be a positive integer!")
  }
  k <- floor(n / l)
  
  ## set seed
  if(!(is.null(seed) || is.na(seed))) set.seed(seed)
  ## bootstrap samples
  res <- replicate(B, 
  {
    j <- sample(1:(n-l+1), k, replace = TRUE)
    x_star <- data[as.vector(sapply(j, function(j) j:(j+l-1)))]
    estim(x_star, version, constant)#, alpha)
  })

  return(list(mean(res > stat), quantile(res, 1 - level)))
}

