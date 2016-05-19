\name{pathCoeff}
\alias{pathCoeff}
\title{Compute path coefficients at a given time lag}
\description{Compute path coefficients.}
\usage{pathCoeff(x, lag = NULL, cumul = FALSE, conf = 0.95)}
\arguments{
  \item{x}{An object of class \code{dlsem}.}
  \item{lag}{A non-negative integer or a vector of non-negative integers indicating the time lag at which path coefficients must be computed.
  If \code{NULL}, all relevant time lags are considered.}
  \item{cumul}{Logical. If \code{FALSE} (the default), instantaneous path coefficients are returned, otherwise cumulative ones are provided.}
  \item{conf}{The confidence level for edges in the object of class \code{dlsem}. Default is 0.95.}
}
\value{A list of matrices, one for each dimension of argument \code{lag}. Each matrix contains point estimate and \code{conf}*100\% confidence interval for each path coefficient.}
%\references{}
\author{Alessandro Magrini <magrini@disia.unifi.it>}
\seealso{\link{dlsem}}
\examples{
data(agres)
mycode <- list(
  GVA~quec(NPATENT,0,4),
  PPI~quec(NPATENT,0,6)+quec(GVA,0,10),
  ENTR_INCOME~quec(NPATENT,0,3)+quec(GVA,1,10)
  )
mod0 <- dlsem(mycode,group="COUNTRY",context=c("GDP","FARM_SIZE"),
  data=agres,uniroot.check=TRUE,imputation=FALSE,log=TRUE)
pathCoeff(mod0,lag=c(0,5,10,15),cumul=TRUE)
}
