estim <- function(x, version, alpha = 0.5)
{
  n <- length(x)
  
  if(version == "empVar")
  {
    res <- .Call("CUSUM", as.numeric(x^2))
  } else if(version == "MD")
  {
    y <- cumstats::cummedian(x)
    res <- .Call("MD", as.numeric(x), as.numeric(y), as.numeric(n)) / (1:(n-1))
  } else if(version == "GMD")
  {
    res <- .Call("GMD", as.numeric(x), as.numeric(n)) / ((1:(n-1)) * (2:n)) * 2
  } else if(version == "MAD")
  {
    res <- sapply(2:n, function(k) mad(x[1:k], constant = constant))
  } else if(version == "Qalpha")
  {
    res <- Qalpha(x, alpha)
  }
  return(max(abs(res - res[n-1]) * 2:n))
}



##'@name dbb 
##'@description Dependent block bootstrap (overlapping, non-circular).
##'
##'@param x (time series)
##'@param l (block length, 1 <= l)
##'@param B (number of bootstrap samples, numeric > 0)
##'@param seed (start for random number generator)
##'       
##'@return 
dbb <- function(stat, data, version, control = list(), alpha = 0.5)
{
  n <- length(x)
  
  browser()
  
  ## argument check
  # if(!is.na(l) && (!is.numeric(l) || l < 1 || l > n))
  # {
  #   stop("l must be a positive integer and cannot be larger than the length of x!")
  # }
  # if(missing(l) | is.na(l)) 
  # {
  #   rho <- abs(cor(x[1:(n-1)], x[2:n], method = "spearman"))
  #   l <- max(ceiling(n^(1/3) * ((2 * rho) / (1 - rho^2))^(2/3)), 1)
  #   l <- tryCatch(as.integer(l), error = function(e) stop("Integer overflow in default l estimation. Please specify a value manually."), 
  #                 warning = function(w) stop("Integer overflow in default l estimation. Please specify a value manually."))
  # }
  
  B <- control$B
  l <- control$l
  seed <- control$seed
  
  if(!is.numeric(B) || B < 1)
  {
    stop("B has to be a positive integer!")
  }
  
  k <- floor(n / l)
  
  ## set seed
  if(!(is.null(seed) || is.na(seed))) set.seed(seed)
  ## bootstrap samples
  res <- replicate(B, 
  {
    j <- sample(1:(n-l+1), k, replace = TRUE)
    x_star <- data[as.vector(sapply(j, function(j) j:(j+l-1)))]
    estim(x_star, version, alpha)
  })
  
  return(mean(res > stat))
}

