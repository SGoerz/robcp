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
  if(!(class(x) %in% c("numeric", "integer") & 
       class(y) %in% c("integer", "numeric"))) 
  {
    stop("x any y have to be numeric!")
  }
  
  if(length(x) == 1 | length(y) == 1)
  {
    return(median(x - y))
  }
  
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