kthPair <- function(x, y, k)
{
  n <- length(x)
  m <- length(y)
  stopifnot(1 <= k && k <= n * m)
  
  x <- sort(x, decreasing = TRUE)
  y <- sort(y, decreasing = TRUE)
  
  res <- .Call("KthPair", as.numeric(x), as.numeric(y), as.numeric(n),
               as.numeric(m), as.numeric(k))
  
  return(res)
}