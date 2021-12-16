## Hodges-Lehmann

##'@name HodgesLehmann
##'@title Hodges Lehmann Test Statistic
##'@description Computes the test statistic for the Hodges-Lehmann change point test.
##'@param x time series (numeric or ts vector).
##'@param b bandwidth for u_hat() (numeric > 0).
##'@param method method for estimating the long run variance.
##'@param control a list of control parameters.
##'@return Test statistic (numeric value) with the attribute cp-location 
##'        indicating at which index a change point is most likely. Is an S3 
##'        object of the class cpStat        
HodgesLehmann <- function(x, b_u = "nrd0", method = "subsampling", control = list(), 
                          p1, p2)
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
  if(length(x) < 2)
  {
    stop("x must consist of at least 2 observations!")
  }
  if(is.null(control$distr) || is.na(control$distr)) 
  {
    control$distr <- TRUE
  }
  if(is.null(control$overlapping)) 
  {
    control$overlapping <- TRUE
  }
  method <- match.arg(method, c("subsampling", "kernel", "bootstrap"))
  ## end argument check
  n <- length(x)
  
  # ## first iteration (k == 1)
  # medDiff <- medianDiff(x[2:n], x[1])
  # x.adj <- x - c(0, rep(medDiff, n - 1))
  # 
  # ## determine b_u adaptively
  # if(is.na(b_u))
  # {
  #   diffs <- unlist(sapply(1:(n-1), function(i)
  #   {
  #     x.adj[(i+1):n] - x.adj[i]
  #   }))
  # 
  #   b_u <- bw.SJ(diffs)
  #   
  #   # diffs <- unlist(sapply(1:(n-1), function(i)
  #   # {
  #   #   temp <- rep(x.adj[(i+1):n], each = i)
  #   #   x.adj[1:i] - temp
  #   # }))
  #   # 
  #   # # sometimes diffs is too large
  #   # b_u <- tryCatch(bw.SJ(diffs), error = function(e) bw.nrd0(diffs))
  # }
  # 
  # ## first Mn
  # Mn <- u_hat(x.adj, b_u, "QS") * (n-1) / n^2 * abs(medDiff)
  
  Mn <- sapply(1:(n-1), function(k)
  {
    medDiff <- medianDiff(x[(k+1):n], x[1:k])
    x.adj <- x - c(rep(0, k), rep(medDiff, n - k))
    #x.adj <- x - c(rep(median(x[1:k]), k), rep(median(x[(k+1):n]), n - k))
    
    diffs <- rep(x.adj, each = n) - as.numeric(x.adj)
    diffs[which(diffs == 0)] = NA
    
    #diffs <- rep(x.adj[1:k], each = n - k) - as.numeric(x.adj[(k+1):n])

    #dens <- u_hat(x.adj, b_u, "QS") 
    dens <- density(diffs, na.rm = TRUE, from = 0, to = 0, n = 1, 
                    bw = b_u)$y
    
    dens * k / n * (1 - k / n) * abs(medDiff)
  })
  
  k <- which.max(Mn)

  if((method == "subsampling" & (is.null(control$l) || is.na(control$l))) | 
     (method == "kernel" & (is.null(control$b_n) || is.na(control$b_n))) | 
     (method == "bootstrap" & (is.null(control$l) || is.na(control$l))))
  {
    n <- length(x)
    x.adj <- x
    #x.adj[(k+1):n] <- x.adj[(k+1):n] - mean(x[(k+1):n]) + mean(x[1:k])
    rho <- abs(cor(x.adj[-n], x.adj[-1], method = "spearman"))
    
    param <- max(ceiling(n^(p1) * ((2 * rho) / (1 - rho^2))^(p2)), 1)
    control$b_n <- min(param, n-1)
    control$l <- control$b_n
  }

  Tn <- sqrt(n) * max(Mn) / sqrt(lrv(x, method = method, control = control))

  attr(Tn, "cp-location") <- k
  class(Tn) <- "cpStat"

  return(Tn)
}

u_hat <- function(x, b_u = "SJ")
{
  n <- length(x)
  diffs <- rep(x, each = n) - as.numeric(x)
  diffs[which(diffs == 0)] = NA
  
  print(b_u)
  
  return(density(diffs, na.rm = TRUE, from = 0, to = 0, n = 1, bw = b_u)$y)
}

# ## outdated?
# u_hat <- function(x, b_u, kFun = "QS")
# {
#   if(b_u <= 0) 
#     stop("b must be numeric, greater than 0 and smaller than the length of the time series!")
#   
#   n <- length(x)
#   kFun <- pmatch(kFun, c("bartlett", "FT", "parzen", "QS", "TH", "truncated",
#                          "Gaussian"))
#   res <- .Call("u_hat", as.numeric(x), as.numeric(b_u), as.numeric(kFun))
#   return(res)
# }