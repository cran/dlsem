\name{edgeCoeff}
\alias{edgeCoeff}
\title{Coefficients associated to each edge}
\description{Compute the coefficient associated to each edge at different time lags.}
\usage{edgeCoeff(x, lag = NULL, conf = 0.95)}
\arguments{
  \item{x}{An object of class \code{dlsem}.}
  \item{lag}{A non-negative integer or a vector of non-negative integers indicating the time lag at which path coefficients must be computed.
  If \code{NULL}, all relevant time lags are considered.}
  \item{conf}{The confidence level for edges in the object of class \code{dlsem}. Default is 0.95.}
}
\value{A list of matrices, one for each dimension of argument \code{lag}. Each matrix contains point estimate and \code{conf}*100\% confidence interval for each path coefficient.}
%\note{}
%\references{}
\author{Alessandro Magrini <magrini@disia.unifi.it>}
\seealso{\link{dlsem}}
\examples{
data(agres)
mycode <- list(
  GVA~quec(NPATENT,1,15),
  PPI~quec(NPATENT,0,13)+quec(GVA,0,14),
  ENTR_INCOME~quec(NPATENT,0,14)+quec(GVA,1,14)
  )
myfit <- dlsem(mycode,group="COUNTRY",context=c("GDP","FARM_SIZE"),
  data=agres,uniroot.check=TRUE,imputation=FALSE,log=TRUE)
edgeCoeff(myfit,lag=0:10)
}
