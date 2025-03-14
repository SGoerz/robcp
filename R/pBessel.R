# pBessel: asymptotic cumulative distribution of the maximum of a h dim. 
#         squared Bessel bridge
pBessel <- function(tn, p) 
{
  ## argument check
  if(any(length(tn) != 1, length(p) != 1, !is.numeric(tn), !is.numeric(p)))
  {
    stop("Incorrect arguments provided!")
  }
  if(p <= 0)
  {
    stop("h has to be a positive integer!")
  }
  tn <- as.numeric(tn)
  ## end argument check
  
  if (p == 1) {
    return(pKSdist(sqrt(tn)))
  }
  if (tn <= 0) return(0)
  if (tn >= (p / 3.5 + 300 / log(p))) return(1)
  
  data("zeros", envir = environment()) 
  tn <- sqrt(tn)
  vfak <- 4 / tn^2
  zerosv <- get("zeros")[, (p - 2) + 1]
  z <- outer(zerosv, 1 / tn)
  fuval <- exp( (p - 2) * log(z) - 1/2 * z^2 - 
                  lgamma(p / 2) - p / 2 * log(2)) / 
    matrix(besselJ(zerosv, p / 2)^2, 
           nrow = 50)
  
  fuval <- apply(fuval, 2, sum)
  return(vfak * fuval)
}