# kthPair <- function(x, y, k)
# {
#   n <- length(x)
#   m <- length(y)
#   stopifnot(1 <= k && k <= n * m)
#   
#   x <- sort(x, decreasing = TRUE)
#   y <- sort(y, decreasing = TRUE)
#   
#   res <- .Call("KthPair", as.numeric(x), as.numeric(y), as.numeric(n),
#                as.numeric(m), as.numeric(k))
#   
#   return(res)
# }

##' Computes the median of the set X - Y. X - Y denotes the set {x - y | x \in X, y \in Y}
##' 
##' @param x,y Numeric vectors
##' 
##' @return The median element of X - Y
##' 
##' @examples 
##' x <- rnorm(100)
##' y <- runif(100)
##' medianDiff(x, y)
medianDiff <- function(x, y)
{
  ## checks missing
  
  nm <- length(x) * length(y)
  
  if(nm %% 2 == 0)
  {
    k <- nm / 2
    return((kthPair(x, -y, k) + kthPair(x, -y, k + 1)) / 2)
  } else
  {
    k <- ceiling(nm / 2)
    return(kthPair(x, -y, k))
  }
}