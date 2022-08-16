##'Print method for change point statistics
##'
##'@description
##'Prints the value of the test statistic and add the most likely change point location.
##'
##'@param x object of the class 'cpStat'.
##'@param ... further arguments passed to or from other methods.
##'
##'@seealso \code{\link{print}}, \code{\link{wmw_stat}}, \code{\link{CUSUM}}, \code{\link{HodgesLehmann}}, 
print.cpStat <- function(x, ...)
{
  loc <- attr(x, "cp-location")
  print(round(as.numeric(x), digits = getOption("digits")), ...)
  cat("location: ", loc)
  return(invisible(x))
}


##'Plot method for change point statistics
##'
##'@description 
##'Plots the trajectory of the test statistic process with a red line indicating 
##'the critical value (default: to the level alpha = 0.05) and a blue line indicating the most
##'probable change point location
##'
##'@param x object of the class 'cpStat'.
##'@param ylim the y limits of the plot.
##'@param xaxt a character which specifies the x axis type (see \code{\link{par}}).
##'@param crit.val critical value of the test.
##'@param ... other graphical parameters (see \code{\link{par}}).
##'
##'@details
##' * Default \code{ylim} is \code{c(min(c(data, 1.358)), max(c(data, 1.358)))}.
##' * Default \code{xaxt} is the similar to the option \code{"s"}, only that there 
##'   is a red labelled tick at the most probable change point location. Ticks too 
##'   close to this will be suppressed.
##'   
##'@seealso \code{\link{plot}}, \code{\link{par}}, \code{\link{CUSUM}}, \code{\link{HodgesLehmann}}, 
##'         \code{\link{wmw_stat}}
plot.cpStat <- function(x, ylim, xaxt, crit.val, ...)
{
  if(missing(crit.val))
  {
    if(!is.null(attr(x, "m"))) 
    {
      m <- attr(x, "m")
      crit.val <- uniroot(Vectorize(function(x) pBessel(x, m) - 0.95), 
                          interval = c(0, m + 1))$root
    }
    else crit.val <- 1.358
  }
  
  data <- attr(x, "teststat")
  if(missing(ylim)) ylim <- c(min(c(data, crit.val)), max(c(data, crit.val)))
  if(missing(xaxt))
  {
    xaxis <- TRUE
    xaxt <- "n"
  }
  plot(data, ylim = ylim, xaxt = "n", ...)
  k <- attr(x, "cp-location")
  abline(v = k, col = "blue")
  
  if(xaxis)
  {
    ticks <- axTicks(1)
    ticks <- ticks[abs(ticks - k) / length(data) >= 0.05]
    axis(1, ticks)
    axis(1, k, col.axis = "blue", col.ticks = "blue", tick = TRUE)
  }
  abline(h = crit.val, col = "red")
  
  return(invisible(NULL))
}