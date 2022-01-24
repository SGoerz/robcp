##'CUSUM: computes the test statistic for the CUSUM test
##'
##'input: y (time series; vector, matrix or ts object)
##'       fun (psi function; character string)
##'       b_n (bandwidth for the long run variance estimation; numeric > 0)
##'       inverse (character string specifying the method of inversion)
##'       ... (further arguments for the inverse-computing functions)
##'       
##'output: test statistic (numeric) with the attribute "cp-location" indicating 
##'        at which index a change point is most likely
##'        -> class "cpStat"

CUSUM <- function(x, method = "kernel", control = list(), inverse = "Cholesky", p1, p2, ...)
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
  method <- match.arg(method, c("subsampling", "kernel", "bootstrap"))
  ## end argument check
  
  
  if(is(x, "matrix"))
  {
    sigma <- lrv(x, method = method, control = control)
    
    m <- ncol(x)
    n <- nrow(x)
    
    if(inverse == "Cholesky")
    {
      # Modified Cholesky Factorization
      mchol <- modifChol(sigma, ...)
      swaps <- attr(mchol, "swaps")
      # invert sigma
      sigma.inv <- chol2inv(mchol) 
    } else if(inverse == "generalized")
    {
      if(!requireNamespace("MASS", quietly = TRUE)) 
      {
        stop("Package \"MASS\" needed for this function to work. Please install it.",
             call. = FALSE)
      }
      sigma.inv <- MASS::ginv(sigma, ...)
      swaps <- 0:(m - 1)
    }
    
    temp <- .Call("CUSUM_ma", as.numeric(x), as.numeric(sigma.inv), 
                 as.numeric(swaps), as.numeric(n), as.numeric(m))
  } else
  {
    temp <- .Call("CUSUM", as.numeric(x))
    k <- which.max(temp)
    
    if((method == "subsampling" & (is.null(control$l) || is.na(control$l))) | 
       (method == "kernel" & (is.null(control$b_n) || is.na(control$b_n))) | 
       (method == "bootstrap" & (is.null(control$l) || is.na(control$l))))
    {
      n <- length(x)
      x.adj <- x
      x.adj[(k+1):n] <- x.adj[(k+1):n] - mean(x[(k+1):n]) + mean(x[1:k])
      rho <- cor(x.adj[-n], x.adj[-1], method = "spearman")
      
      param <- max(ceiling(n^(p1) * ((2 * rho) / (1 - rho^2))^(p2)), 1)
      param <- min(param, n-1)

      control$b_n <- param
      control$l <- param
    }
    
    if(method == "kernel" & (is.null(control$kFun) || is.na(control$kFun)))
    {
      control$kFun <- "TH"
    }
    
    sigma <- sqrt(lrv(x, method = method, control = control))
    temp <- temp / sigma
  }
  
  erg <- max(temp)
  attr(erg, "cp-location") <- k
  attr(erg, "data") <- ts(x)
  attr(erg, "lrv-estimation") <- method
  attr(erg, "sigma") <- sigma
  if(method == "kernel") attr(erg, "b_n") <- control$b_n else
    attr(erg, "l") <- control$l
  attr(erg, "teststat") <- ts(temp)
  
  class(erg) <- "cpStat"
  
  return(erg)
}

##'print.cpStat: print method for change point statistics
##'              prints the value of the test statistic and add the most likely
##'              change point location
##'@name CUSUM
print.cpStat <- function(x, ...)
{
  loc <- attr(x, "cp-location")
  print(round(as.numeric(x), digits = getOption("digits")), ...)
  cat("location: ", loc)
  return(invisible(x))
}


plot.cpStat <- function(x, ylim, xaxt, ...)
{
  data <- attr(x, "teststat")
  if(missing(ylim)) ylim <- c(min(c(data, 1.358)), max(c(data, 1.358)))
  if(missing(xaxt)) xaxt <- "n"
  plot(data, ylim = ylim, xaxt = "n", ...)
  k <- attr(x, "cp-location")
  abline(v = k, col = "blue")
  
  ticks <- axTicks(1)
  ticks <- ticks[abs(ticks - k) / length(data) >= 0.05]
  axis(1, ticks)
  axis(1, k, col.axis = "blue", col.ticks = "blue", tick = TRUE)
  abline(h = 1.358, col = "red")
  return(invisible(NULL))
}
