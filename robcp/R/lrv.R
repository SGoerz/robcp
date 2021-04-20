##'lrv: estimates the long run variance resp. covariance matrix of the supplied
##'     time series
##'
##'input: x (time series; numeric vector, matrix or ts object)
##'       b_n (bandwidth for kernel-based estimation; numeric; 
##'            default: length of time series ^ (1/3))
##'       l (block length for subsampling or bootstrap; numeric; 
##'          default: ??????)
##'       method (long run variance estimation method [for the univariate case];
##'               one of "kernel", "subsampling", "bootstrap")
##'       B (number of bootstrap samples; numeric)
##'       gamma0 (for the kernel-based estimation: if the estimated lrv is <= 0,
##'               should only the estimated value to the lag 0 be returned?; 
##'               default: TRUE)
##'       
##'output: long run variance (numeric value) or long run covariance matrix 
##'        (numeric matrix with dim. m x m, when m is the number of columns)

lrv <- function(x, b_n, l, method = "kernel", B, gamma0 = TRUE)
{
  ## argument check
  if(is(x, "ts"))
  {
    class(x) <- "numeric"
  } 
  if(!(is(x, "matrix") || is(x, "numeric") || is(x, "integer")))
  {
    stop("x must be a numeric or integer vector or matrix!")
  }
  if(!missing(b_n) && (!is(b_n, "numeric") || b_n <= 0))
  {
    stop("b_n must be numeric and greater than 0!")
  }
  if(!missing(l) && (!is(l, "numeric") || l <= 0))
  {
    stop("l must be numeric and greater than 0!")
  }
  ## end argument check
  
  if(is(x, "matrix"))
  {
    m <- ncol(x)
    n <- nrow(x)

    if(missing(b_n)) b_n <- n^(1/3)
    if(b_n > n)
      stop("The bandwidth b_n cannot be larger than the length of the time series!")
    
    x_cen <- apply(x, 2, function(x) x - mean(x))
    
    erg <- .Call("lrv_matrix", as.numeric(x_cen), 
                 as.numeric(n), as.numeric(m), as.numeric(b_n), 
                 PACKAGE = "robcp")
    
    erg <- matrix(erg, ncol = m)
  } else
  {
    method <- match.arg(method, c("subsampling", "kernel", "bootstrap"))
    erg <- switch(method, 
           "kernel" = lrv_kernel(x, b_n, gamma0), 
           "subsampling" = lrv_subs(x, l), 
           "bootstrap" = lrv_dwb(x, l, B))
  }
  
  return(erg)
}


##'kernel estimation
##'
##'input: x (time series)
##'       b_n (bandwidth)
##'       gamma0 (use only hat(gamma)(0) when estimation <= 0?)
lrv_kernel <- function(x, b_n, gamma0 = TRUE)
{
  n <- length(x)
  if(missing(b_n)) b_n <- n^(1/3)
  if(b_n > n) 
    stop("The bandwidth b_n cannot be larger than the length of the time series!")
  
  x_cen <- x - mean(x)
  erg <- .Call("lrv", as.numeric(x_cen), as.numeric(b_n),
               PACKAGE = "robcp")
  if(erg < 0 & gamma0)
  {
    warning("Estimated long run variance was < 0; only the estimated autocovariance to lag 0 is returned!")
    erg <- (n - 1) / n * var(x)
  }
  return(erg)
}

##'overlappign subsampling estimation
##'
##'input: x (time series)
##'       l (block length; numeric; 1 <= l <= length(x))
lrv_subs <- function(x, l)
{
  ecdf.values <- ecdf(x)(x)
  
  res <- .Call("lrv_subs", as.numeric(ecdf.values), as.numeric(l))
  return(res)
}


##'dependent wild bootstrap estimation 
##'
##'input: x (time series)
##'       l (block length??, 1 <= l)
##'       B (number of bootstrap samples, numeric > 0)
lrv_dwb <- function(x, l, B)
{
  browser()
  n <- length(x)
  
  if(!missing(l) && (!is.numeric(l) || l < 1 || l > n))
  {
    stop("l must be a positive integer and cannot be larger than the length of x!")
  }
  if(missing(l)) l <- sqrt(n) ##????????????????
  if(!missing(B) && (!is.numeric(B) || B < 1)) 
  {
    stop("B has to be a positive integer!")
  }
  if(missing(B)) B <- 1000
  
  sigma.ma <- matrix(.Call("gen_matrix", as.numeric(n), as.numeric(l)), ncol = n)
  
  # sigma.ma <- as.matrix(dist(0:(n-1), diag = TRUE, upper = TRUE)) / l
  # sigma.ma <- apply(sigma.ma, 1, K)
  
  if(!requireNamespace("pracma", quietly = TRUE)) 
  {
    stop("Package \"pracma\" needed for the dependent wild bootstrap to work. Please install it.",
         call. = FALSE)
  }

  ## dependency matrix
  sigma.root <- pracma::sqrtm(sigma.ma)$B
  ## bootstrap samples
  dwb <- replicate(B, {
    z <- rnorm(n)
    eps <- sigma.root %*% z
    x_star <- mean(x) + (x - mean(x)) * eps
    mean(x_star)
  })
  
  return(var((dwb - mean(x)) * sqrt(n)))
}

