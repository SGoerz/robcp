opt.param <-  function(x)
{
  n <- length(x)
  if(n <= 5) stop("For automatic bandwidth selection x must consist of at least 6 observations!")
  kappa <- max(5, sqrt(log10(n)))
  m <- round(n^(1/3), 5)
  rho <- abs(acf(x, plot = FALSE, lag.max = m + kappa)[[1]][, , 1][-1])
  i <- 1
  cond <- 2 * sqrt(log10(n) / n)
  repeat
  {
    if(i > m - 1 || max(rho[i:(kappa+i)]) <= cond) break
    i <- i + 1
  }
  
  return(i)
}

## Tests for Scale Changes Based on Pairwise Differences (Gerstenberger et. al.)

##'@name scale_stat
##'@title Test statistic to detect Scale Changes 
##'@description Computes the test statistic for CUSUM-based tests on scale changes.
##'@param x time series (numeric or ts vector).
##'@param version variance estimation method. One of "empVar", "MD", "GMD", "MAD", "Qalpha".
##'@param control a list of control parameters.
##'@param constant scale factor for the MAD.
##'@param alpha quantile of the distribution function of all absolute pairwise differences used in \code{version = "Qalpha"}.
##'@return Test statistic (numeric value) with the attribute cp-location 
##'        indicating at which index a change point is most likely. Is an S3 
##'        object of the class cpStat   

scale_stat <- function(x, version = "empVar", method = "kernel", control = list(), 
                       constant = 1.4826, alpha = 0.5)
{
  ## argument check
  if(is(x, "ts"))
  {
    class(x) <- "numeric"
    tsp <- attr(x, "tsp")
  }
  if(!(is(x, "matrix") || is(x, "numeric") || is(x, "integer")))
  {
    stop("x must be a numeric or integer vector or matrix!")
  }
  n <- length(x)
  if(n < 2)
  {
    stop("x must consist of at least 2 observations!")
  }
  if(is.null(control$kFun) || is.na(control$kFun))
  {
    control$kFun <- "SFT"
  }
  
  if(method == "kernel" && (is.null(control$b_n) || is.na(control$b_n)))
  {
    control$b_n <- max(opt.param(x), opt.param(x^2))
  }
  
  ## end argument check
  
  control$version <- version
  
  if(version == "empVar")
  {
    control$loc <- mean(x)
    control$scale <- var(x)

    stat <- .Call("CUSUM", as.numeric(x^2))
    # sigma <-  sqrt(lrv(x, method = "kernel", control = control))
    res <- max(stat)
    k <- which.max(res)
  } else 
  {
    if(version == "MD")
    {
      if(!requireNamespace("cumstats", quietly = TRUE)) 
      {
        stop("Package \"cumstats\" needed for the MD. Please install it.",
             call. = FALSE)
      }
      y <- cumstats::cummedian(x)
      res <- .Call("MD", as.numeric(x), as.numeric(y), as.numeric(n)) / (1:(n-1))

      control$loc <- median(x)
    } else if(version == "GMD")
    {
      res <- .Call("GMD", as.numeric(x), as.numeric(n)) / ((1:(n-1)) * (2:n)) * 2
    } else if(version == "MAD")
    {
      res <- sapply(2:n, function(k) mad(x[1:k], constant = constant))
      control$loc <- median(x)
    } else if(version == "Qalpha")
    {
      res <- Qalpha(x, alpha)
      # sorted <- x[1]
      # res <- sapply(2:n, function(k)
      # {
      #   sorted <<- sort(c(sorted, x[k]), decreasing = TRUE)
      #   a <- ceiling(k * (k - 1) / 2 * (1 - beta))
      #   # a <- floor(k * (k - 1) / 2 * (1 - beta)) + 1
      #   kthPair(sorted[1:(k-1)], -sorted[2:k], a)
      # })
      
      control$loc <- alpha
    } else 
    {
      stop("version not supported.")
    }

    control$scale <- res[n-1]
    control$distr <- FALSE
    
    stat <- res - res[n-1]
    stat <- 2:n * abs(stat)
    k <- which.max(stat) + 1
    res <- max(stat) / sqrt(n)
  } 
  
  if(method == "kernel") 
  {
    sigma <- sqrt(lrv(x, method = "kernel", control = control))
    res <- res / sigma 
    
    attr(res, "lrv-estimation") <- "kernel"
    attr(res, "sigma") <- sigma
    attr(res, "param") <- control$b_n 
    attr(res, "kFun") <- control$kFun
    stat <- stat / sigma
    
  } else if(method == "bootstrap")
  {
    attr(res, "lrv-estimation") <- NA
  } 
  
  attr(res, "cp-location") <- as.integer(k)
  attr(res, "teststat") <- ts(stat)
  class(res) <- "cpStat"
  
  return(res)
}



## Scale change test

##'@name scale_cusum
##'@title Tests for Scale Changes Based on Pairwise Differences
##'@description Performs the CUSUM-based test on changes in the scale.
##'@param x time series (numeric or ts vector).
##'@param version variance estimation method. One of "empVar", "MD", "GMD", "MAD", "QBeta".
##'@param control a list of control parameters.
##'@param constant scale factor for the MAD.
##'@param beta quantile of the distribution function of all absolute pairwise differences used in \code{version = "QBeta"}.
##'@param fpc finite population correction (boolean).
##'@param tol tolerance of the distribution function (numeric), which is used do compute p-values.
##'@return A list fo the class "htest" containing

scale_cusum <- function(x, version = "empVar", method = "kernel", control = list(), 
                        constant = 1.4826, alpha = 0.5, fpc = TRUE, tol = 1e-8)
{
  Dataname <- deparse(substitute(x))
  
  stat <- scale_stat(x = x, version = version, method = method, control = control,
                     constant = constant, alpha = alpha)
  location <- attr(stat, "cp-location")
  names(stat) <- "S"
  
  if(fpc) stat <- stat + 1.46035 / sqrt(2 * pi) / sqrt(length(x))
  
  if(method == "kernel")
  {
    p.val <- 1 - pKSdist(stat, tol)
    erg2 <- list(lrv = list(method = "kernel", 
                            param = attr(stat, "param"), 
                            value = attr(stat, "sigma")))
  } else if(method == "bootstrap")
  {
    control$l <- opt.param(x)
    if(is.null(control$B)) control$B <- 1 / tol
    p.val <- dbb(stat, data = x, version = version, control = control, alpha = alpha)
    erg2 <- list(bootstrap = list(param = control$l))
  } else
  {
    stop("method not supported.")
  }
  
  erg <- list(alternative = "two-sided", method = "CUSUM test for scale changes",
              data.name = Dataname, statistic = stat,
              p.value = p.val,
              cp.location = location)
  
  erg <- c(erg, erg2)
  
  class(erg) <- "htest"
  return(erg)
}

