\name{plot.cpStat}
\alias{plot.cpStat}
\title{Plot method for change point statistics}
\usage{
\method{plot}{cpStat}(x, ylim, xaxt, crit.val, ...)
}
\arguments{
\item{x}{object of the class 'cpStat'.}

\item{ylim}{the y limits of the plot.}

\item{xaxt}{a character which specifies the x axis type (see \code{\link{par}}).}

\item{crit.val}{critical value of the test. Default: 1.358.}

\item{...}{other graphical parameters (see \code{\link{par}}).}
}
\description{
Plots the trajectory of the test statistic process together with a red line indicating 
the critical value (alpha = 0.05) and a blue line indicating the most
probable change point location
}
\details{
\itemize{
\item Default for \code{ylim} is \code{c(min(c(data, 1.358)), max(c(data, 1.358)))}.
\item Default for \code{xaxt} is the simliar to the option \code{"s"}, only that there 
  is a red labelled tick at the most probable change point location. Ticks too 
  close to this will be suppressed.
}
}
\seealso{
\code{\link{plot}}, \code{\link{par}}, \code{\link{CUSUM}}, \code{\link{HodgesLehmann}}, 
        \code{\link{wilcox_stat}}
}


\author{
Sheila \enc{Görz}{Goerz}
}
