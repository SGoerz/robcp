# not in use at the moment
# fastMedian <- function(x)
# {
#   if(any(is.na(x))) return(NA)
#   erg <- .Call("fastMedian", as.numeric(x))
#   return(erg)
# }
# 
# fastMAD <- function(x, med = NA, constant = 1.4826)
# {
#   if(any(is.na(x))) return(NA)
#   if(is.na(med)){med <- fastMedian(x)}
#   
#   erg <- .Call("fastMAD", as.numeric(x), as.numeric(med), as.numeric(constant))
#   return(erg)
# }

# argCheck <- function(y)
# {
#   if(is(y, "ts"))
#   {
#     timeseries <- TRUE
#     tsp <- attr(y, "tsp")
#   } else
#   {
#     timeseries <- FALSE
#     tsp <- NULL
#   }
#   if(!(is(y, "matrix") | is(y, "numeric") | is(y, "integer")))
#   {
#     stop("x must be numeric and either be a vector, a matrix or a ts object!")
#   }
#   
#   return(list(timeseries = timeseries, tsp = tsp))
# }


# # HLm: marginal huberized location
# # HLg: global huberized location
# # SLm: marginal sign location
# # SLg: global sign location
# # HCm: marginal huberized covariance
# # HCg: global huberized covariance
# # SCm: marginal sign covariance
# # SCg: global sign covariance
# psi <- function(y, fun = "HLm", k, constant = 1.4826)
# {
#   if(!(is(y, "matrix") | is(y, "numeric") | is(y, "integer")))
#   {
#     stop("x must be numeric and either be a vector or a matrix!")}
# 
#   fun <- which(c("HLm", "HLg", "SLm", "SLg", "HCm", "HCg", "SCm", "SCg") == fun)
#   if(length(fun) == 0) stop("Wrong argument to \'fun\'")
# 
#   if(missing(k))
#   {
#     if(fun == 2 || fun == 6)
#     {
#       k <- sqrt(qchisq(0.8, df = ifelse(is(y, "matrix"), ncol(y), 1)))
#     }
#     else k <- 1.5
#   }
#   
#   if(is(y, "matrix"))
#   {
#     m <- ncol(y)
#     n <- nrow(y)
#     med <- apply(y, 2, median)
#     MAD <- apply(y, 2, function(x) mad(x, constant = constant))
#   }
#   else
#   {
#     if(fun %in% 7:8) 
#     {
#       stop("VCm and VCg are not available for a one-dimensional timeline")
#     }
#     n <- length(y)
#     m <- 1
#     med <- median(y)
#     MAD <- mad(y, constant = constant)
#   }
#   
#   erg <- .Call("psi", as.numeric(y), as.numeric(fun), as.numeric(n), 
#                as.numeric(m), as.numeric(k), #as.numeric(constant), 
#                # median and mad are being computed in R 
#                as.numeric(med), as.numeric(MAD))
#   
#   if(is(y, "matrix")) erg <- matrix(erg, nrow = n)
#   
#   return(erg)
# }
# 
# 
# h_cumsum <- function(y, fun = "HLm", k, constant)
# {
#   if(!(is(y, "matrix") || is(y, "numeric") || is(y, "integer")))
#   {
#     stop("x must be a numeric or integer vector or matrix!")
#   }
# 
#   x <- psi(y, fun, k, constant)
# 
#   if(is(x, "matrix"))
#   {
#     m <- ncol(x)
#     n <- nrow(x)
# 
#     erg <- .Call("cumsum_ma", as.numeric(x), as.numeric(n), as.numeric(m))
#     erg <- matrix(erg, ncol = m)
#   }
#   else
#   {
#     erg <- .Call("h_cumsum", as.numeric(x))
#   }
# 
#   return(erg)
# }


# sigma2 <- function(x, b_n)
# {
#   if(!(is(x, "matrix") || is(x, "numeric") || is(x, "integer")))
#   {
#     stop("x must be a numeric or integer vector or matrix!")
#   }
#   if(!missing(b_n) && b_n <= 0) stop("b_n must be greater than 0")
#   
#   if(is(x, "matrix"))
#   {
#     m <- ncol(x)
#     n <- nrow(x)
#     
#     if(missing(b_n)) b_n <- n^(1/3)
#     
#     x_cen <- apply(x, 2, function(x) x - mean(x))
#     
#     erg <- .Call("sigma_matrix", as.numeric(x_cen), 
#                  as.numeric(n), as.numeric(m), as.numeric(b_n))
#     
#     return(matrix(erg, ncol = m))
#   }
#   
#   if(missing(b_n)) b_n <- length(x)^(1/3)
#   x_cen <- x - mean(x)
#   erg <- .Call("sigma2", as.numeric(x_cen), as.numeric(b_n),
#                PACKAGE = "robcp")
#   return(erg)
# }


# modifChol <- function(x, tau = .Machine$double.eps^(1/3),
#                        tau_bar = .Machine$double.eps^(2/3), mu = 0.1)
# {
#   if(!is(x, "matrix")) stop("x must be a matrix!")
#   
#   n <- nrow(x)
#   if(n != nrow(x)) stop("x must be a square matrix!")
#   
#   erg <- matrix(.Call("cholesky", as.numeric(x), as.numeric(n),as.numeric(tau),
#                as.numeric(tau_bar), as.numeric(mu)), nrow = n)
#   
#   swaps <- erg[, n + 1]
#   erg <- erg[, -(n + 1)]
#   attr(erg, "swaps") <- swaps
#   return(erg)
# }


# teststat <- function(y, fun = "HLm", b_n, k, constant)
# {
#   if(!(is(y, "matrix") || is(y, "numeric") || is(y, "integer")))
#   {
#     stop("x must be a numeric or integer vector or matrix!")
#   }
#   
#   x <- psi(y, fun, k, constant)
#   
#   if(is(x, "matrix"))
#   {
#     sigma <- lrv(x, b_n)
#     
#     # Modified Cholesky Factorization
#     mchol <- modifChol(sigma)
#     swaps <- attr(mchol, "swaps")
#     # invert sigma
#     sigma <- chol2inv(mchol)  
#     
#     m <- ncol(x)
#     n <- nrow(x)
#     
#     erg <- .Call("teststat_ma", as.numeric(x), as.numeric(sigma), 
#                  as.numeric(swaps), as.numeric(n), as.numeric(m))
#     
#     return(erg)
#   }
#    
#   sigma <- sqrt(lrv(x, b_n))
#   erg <- .Call("teststat", as.numeric(x), as.numeric(sigma))
#   return(erg)
# }

# pKSdist <- function(tn, tol = 1e-8)
# {
#   erg <- .Call("pKSdist", as.numeric(tn), as.numeric(tol))
#   return(erg)
# }

# pbessel3 <- function(tn, h) 
# {
#   if (h == 1) {
#     return(pKSdist(sqrt(tn)))
#   }
#   if (tn <= 0) return(0)
#   if (tn >= (h / 3.5 + 300 / log(h))) return(1)
#   
#   data("zeros", envir = environment()) 
#   tn <- sqrt(tn)
#   vfak <- 4 / tn^2
#   zerosv <- get("zeros")[,(h - 2) + 1]
#   z <- outer(zerosv, 1 / tn)
#   fuval <- exp( (h - 2) * log(z) - 1/2 * z^2 - 
#                   lgamma(h / 2) - h / 2 * log(2)) / 
#     matrix(besselJ(zerosv, h / 2)^2, 
#            nrow = 50)
#   
#   fuval <- apply(fuval, 2, sum)
#   return(vfak * fuval)
# }

# huber_cusum <- function(x, fun = "HLm", tol = 1e-8, b_n, k, constant)
# {
#   if(!(is(x, "matrix") || is(x, "numeric") || is(x, "integer")))
#   {
#     stop("x must be a numeric or integer vector or matrix!")
#   }
#   
#   Dataname <- deparse(substitute(x))
#   stat <- teststat(x, fun, b_n, k, constant)
#   names(stat) <- "S"
#   
#   if(is(x, "matrix"))
#   {
#     h <- ncol(x)
#     h <- switch(fun, HCm = h * (h + 1) / 2, HCg = h * (h + 1) / 2, 
#                 VCm = h * (h - 1) / 2, VCg = h * (h + 1) / 2 - 1, h)
#     
#     stat <- (sqrt(stat) + 1.46035 / sqrt(2 * pi) / sqrt(nrow(x)))^2
#     
#     erg <- list(alternative = "two-sided", method = "Huberized CUSUM test",
#                 data.name = Dataname, statistic = stat,
#                 p.value = 1 - pbessel3(stat, h))
#   }
#   else
#   {
#     stat <- stat + 1.46035 / sqrt(2 * pi) / sqrt(length(x))
#     
#     erg <- list(alternative = "two-sided", method = "Huberized CUSUM test",
#                 data.name = Dataname, statistic = stat,
#                 p.value = 1 - pKSdist(stat, tol))
#   }
#   
#   class(erg) <- "htest"
#   return(erg)
# }

