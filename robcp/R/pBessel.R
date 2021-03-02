##'pBessel: asymptotic cumulative distribution of the maximum of a h dim. 
##'         squared Bessel bridge
##'
##'input: tn (test statistic; numeric)
##'       h (dimension of the time series; integer)
##'       
##'output: P(t_n(X) <= tn) (numeric)

pbessel3 <- function(tn, h) 
{
  ## argument check
  if(any(length(tn) != 1, length(h) != 1, !is.numeric(tn), !is.integer(h)))
  {
    stop("Incorrect arguments provided!")
  }
  if(h <= 0)
  {
    stop("h has to be a positive integer!")
  }
  ## end argument check
  
  if (h == 1) {
    return(pKSdist(sqrt(tn)))
  }
  if (tn <= 0) return(0)
  if (tn >= (h / 3.5 + 300 / log(h))) return(1)
  
  data("zeros", envir = environment()) 
  tn <- sqrt(tn)
  vfak <- 4 / tn^2
  zerosv <- get("zeros")[,(h - 2) + 1]
  z <- outer(zerosv, 1 / tn)
  fuval <- exp( (h - 2) * log(z) - 1/2 * z^2 - 
                  lgamma(h / 2) - h / 2 * log(2)) / 
    matrix(besselJ(zerosv, h / 2)^2, 
           nrow = 50)
  
  fuval <- apply(fuval, 2, sum)
  return(vfak * fuval)
}