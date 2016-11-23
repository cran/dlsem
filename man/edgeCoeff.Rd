\name{edgeCoeff}
\alias{edgeCoeff}
\title{Coefficients associated to each edge}
\description{Compute the coefficient associated to each edge at different time lags.}
\usage{edgeCoeff(x, lag = NULL, conf = 0.95)}
\arguments{
  \item{x}{An object of class \code{dlsem}.}
  \item{lag}{A non-negative integer or a vector of non-negative integers indicating the time lag at which coefficients must be computed.
  If \code{NULL}, all relevant time lags are considered.}
  \item{conf}{The confidence level for edges in the object of class \code{dlsem}. Default is 0.95.}
}
\value{A list of matrices, one for each dimension of argument \code{lag}. Each matrix contains point estimate and \code{conf}*100\% confidence interval for each coefficient.}
%\note{}
%\references{}
%\author{Alessandro Magrini <magrini@disia.unifi.it>}
\seealso{\link{dlsem}}
\examples{
data(industry)
mycode <- list(
  Consum~quec(Job,0,6),
  Pollution~quec(Job,1,11)+quec(Consum,1,6)
  )
myfit <- dlsem(mycode,group="Region",exogenous=c("Population","GDP"),
  data=industry,uniroot.check=TRUE,log=TRUE)
edgeCoeff(myfit,lag=0:10)
}
