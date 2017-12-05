\name{edgeCoeff}
\alias{edgeCoeff}
\title{Coefficient associated to each edge at different time lags}
\description{Return point estimation and confidence interval for the coefficient associated to each edge at different time lags.}
\usage{edgeCoeff(x, lag = NULL, conf = 0.95)}
\arguments{
  \item{x}{An object of class \code{dlsem}.}
  \item{lag}{A non-negative integer or a vector of non-negative integers indicating the time lags at which coefficients must be computed.
  If \code{NULL}, all relevant time lags are considered.}
  \item{conf}{The confidence level. Default is 0.95.}
}
\value{A list of matrices, one for each dimension of argument \code{lag}. Each matrix contains point estimate and confidence interval (at level \code{conf}) for each coefficient.}
%\note{}
%\references{}
\seealso{\link{dlsem}.}
\examples{
data(industry)
mycode <- list(
  Consum~quec(Job,0,5),
  Pollution~quec(Job,1,8)+quec(Consum,1,6)
  )
myfit <- dlsem(mycode,group="Region",exogenous=c("Population","GDP"),data=industry,log=TRUE)
edgeCoeff(myfit,lag=0:10)
}
