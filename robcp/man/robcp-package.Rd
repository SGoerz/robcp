\docType{package}
\name{robcp-package}
\alias{robcp}
\alias{robcp-package}
\title{robcp: Robust Change-Point Tests}
\description{
The package includes robust change point tests designed to identify structural changes in time series. 
As external influences might affect the characteristics of a time series, 
such as location, variability, and correlation structure, changes caused by them
and the corresponding time point need to be detected reliably.
Standard methods can struggle when dealing with heavy-tailed data or outliers; 
therefore, this package comprises robust methods that effectively manage extreme
values by either ignoring them or assigning them less weight. 
Examples of these robust techniques include the Median, Huber M-estimator, mean deviation, Gini's mean difference, and \eqn{Q^{\alpha}}.


The package contains the following tests and test statistics:

\strong{Tests on changes in the location}
\describe{
\item{- Huberized CUSUM test: \code{\link{huber_cusum}} (test), \code{\link{CUSUM}} (CUSUM test statistic), \code{\link{psi_cusum}} (transformation function).}{}
\item{- Hodges-Lehmann test: \code{\link{hl_test}} (test), \code{\link{HodgesLehmann}} (test statistic).}{}
\item{- Wilcoxon-Mann-Whitney change point test: \code{\link{wmw_test}} (test), \code{\link{wilcox_stat}} (test statistic).}{}
}

\strong{Tests on changes in the variability}
\describe{
\item{Estimating the variability using the ordinary variance, mean deviation, Gini's mean difference and \eqn{Q^{\alpha}}:}{}
\item{\code{\link{scale_cusum}} (test), \code{\link{scale_stat}} (test statistic).}{}
}

\strong{Tests on changes in the correlation}
\describe{
\item{Estimating the correlation by a sample version of Kendall's tau or Spearman's rho:}{}
\item{\code{\link{cor_cusum}} (test), \code{\link{cor_stat}} (test statistic).}{}
}
}



\author{
\strong{Maintainer}: Sheila \enc{Görz}{Goerz} \email{sheila.goerz@tu-dortmund.de} 

\strong{Author}: Alexander \enc{Dürre}{Duerre} \email{a.m.durre@math.leidenuniv.nl}

\strong{Thesis Advisor}: Roland Fried \email{msnat@statistik.tu-dortmund.de}
}