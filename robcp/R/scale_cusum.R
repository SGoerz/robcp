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
##'@param version variance estimation method. One of "empVar", "MD", "GMD", "Qalpha".
##'@param control a list of control parameters.
##'@param alpha quantile of the distribution function of all absolute pairwise differences used in \code{version = "Qalpha"}.
##'@return Test statistic (numeric value) with the attribute cp-location 
##'        indicating at which index a change point is most likely. Is an S3 
##'        object of the class cpStat   

scale_stat <- function(x, version = c("empVar", "MD", "GMD", "Qalpha"), 
                       method = "kernel", control = list(),# constant = 1.4826,
                       alpha = 0.8)
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
    control$kFun <- "quadratic"
  }
  version <- match.arg(version)
  ## end argument check
  
  control$version <- version
  
  if(version == "empVar")
  {
    stat <- CUSUM_var_cpp(x, x^2)
    # stat <- .Call("CUSUM_var", as.numeric(x), as.numeric(x^2))
    res <- max(stat)
    k <- which.max(stat)
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
      res <- MD_cpp(x, y) / (1:(n-1))
      # res <- .Call("MD", as.numeric(x), as.numeric(y), as.numeric(n)) / (1:(n-1))
    } else if(version == "GMD")
    {
      res <- GMD_cpp(x) / ((1:(n-1)) * (2:n)) * 2
      # res <- .Call("GMD", as.numeric(x), as.numeric(n)) / ((1:(n-1)) * (2:n)) * 2
    # } else if(version == "MAD")
    # {
    #   res <- sapply(2:n, function(k) mad(x[1:k], constant = constant))
    } else if(version == "Qalpha")
    {
      res <- Qalpha(x, alpha)
      control$alpha_Q <- alpha
    } else 
    {
      stop("version not supported.")
    }
    
    control$distr <- FALSE
    
    stat <- res - res[n-1]
    stat <- 2:n * abs(stat) / sqrt(n)
    k <- which.max(stat) + 1
    res <- max(stat) 
  } 
  
  if(method == "kernel" || method == "subsampling") 
  {
    if((method == "kernel" && (is.null(control$b_n) || is.na(control$b_n))) ||
       (method == "subsampling" && (is.null(control$l_n) || is.na(control$l_n))))
    {
      x.adj <- x
      x.adj <- x
      if(k > 1 & k + 1 < n) x.adj[(k+1):n] <- x.adj[(k+1):n] /
        sd(x.adj[(k+1):n]) * sd(x.adj[1:k])
      rho1 <- abs(cor(x.adj[-n], x.adj[-1], method = "spearman"))
      rho2 <- abs(cor((x.adj[-n])^2, (x.adj[-1])^2, method = "spearman"))
      
      ##
      p1 <- 0.5
      p2 <- 0.3
      ##
      
      param <- max(n^(p1) * ((2 * rho1) / (1 - rho1^2))^(p2),
                   n^(p1) * ((2 * rho2) / (1 - rho2^2))^(p2), 1)
      param <- min(param, n-1)
      if(is.na(param)) param <- 1
      control$b_n <- param
      control$l_n <- param
    }
    
    sigma <- sqrt(lrv(x, method = method, control = control))
    res <- res / sigma 
    
    attr(res, "lrv-estimation") <- method
    attr(res, "sigma") <- sigma
    if(method == "kernel")
    {
      attr(res, "kFun") <- control$kFun
      attr(res, "param") <- control$b_n 
    } else
    {
      attr(res, "param") <- control$l_n
    }
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
##'@param version variance estimation method. One of "empVar", "MD", "GMD", "Qalpha".
##'@param control a list of control parameters.
##'@param alpha quantile of the distribution function of all absolute pairwise differences used in \code{version = "Qalpha"}.
##'@param fpc finite population correction (boolean).
##'@param tol tolerance of the distribution function (numeric), which is used do compute p-values.
##'@return A list fo the class "htest" containing
##'
scale_cusum <- function(x, version = c("empVar", "MD", "GMD", "Qalpha"),
                        method = "kernel", control = list(), #constant = 1.4826, 
                        alpha = 0.8, fpc = TRUE, tol, plot = FALSE, level = 0.05)
{
  if(missing(tol))
  {
    if(method == "kernel") tol <- 1e-8 else tol <- 1e-3
  }
  
  Dataname <- deparse(substitute(x))
  
  stat <- scale_stat(x = x, version = version, method = method, control = control,
                     #constant = constant, 
                     alpha = alpha)
  location <- attr(stat, "cp-location")
  names(stat) <- "S"
  
  if(fpc) stat <- stat + 1.46035 / sqrt(2 * pi) / sqrt(length(x))
  
  if(method == "kernel")
  {
    p.val <- 1 - pKSdist(stat, tol)
    erg2 <- list(lrv = list(method = "kernel", 
                            param = attr(stat, "param"), 
                            value = attr(stat, "sigma")))
    if(plot) plot(stat)
  # } else if(method == "bootstrap")
  # {
  #   if(is.null(control$l)) control$l <- opt.param(x)
  #   if(is.null(control$B)) control$B <- 1 / tol
  #   y <- dbb(stat, data = x, version = match.arg(version), control = control,
  #            alpha = alpha, 
  #            #constant = constant,
  #            level = level)
  #   p.val <- y[[1]]
  #   erg2 <- list(bootstrap = list(param = control$l, crit.value = y[[2]]))
  #   if(plot) plot(stat, crit.val = y[[2]])
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

