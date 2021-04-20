## Hodges-Lehmann
## verschiedene Varianzschaetzer

HodgesLehmann <- function(x, b, l)
{
  n <- length(x)
  Mn <- sqrt(n) * max(sapply(1:(n-1), function(k)
  {
    u_hat(x, b) * k / n * (1 - k / n) * abs(medianDiff(x[(k+1):n], x[1:k])) 
  }))
  
  Tn <- Mn / lrv_subs(x, l)
  
  return(Tn)
}


## braucht man das als eigenstaendige Funktion??????
## Verschiedene Kernel-Funktionen?????
u_hat <- function(x, b)
{
  n <- length(x)
  res <- .Call("u_hat", as.numeric(x), as.numeric(b))
  return(res)
}

lrv_subs <- function(x, l)
{
  ecdf.values <- ecdf(x)(x)

  res <- .Call("lrv_subs", as.numeric(ecdf.values), as.numeric(l))
  return(res)
}