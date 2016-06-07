\name{pathAnal}
\alias{pathAnal}
\title{Perform path analysis at different time lags}
\description{Compute the coefficient associated to each path between one or more starting variables and one ending variable at one or more time lags.}
\usage{pathAnal(x, from = NULL, to = NULL, lag = NULL, cumul = FALSE, conf = 0.95, nitt = 50000)}
\arguments{
  \item{x}{An object of class \code{dlsem}.}
  \item{from}{The name of the starting variable, or a vector of names of starting variables.}
  \item{to}{The name of the ending variable.}
  \item{lag}{A non-negative integer or a vector of non-negative integers indicating the time lag at which path coefficients must be computed. If \code{NULL}, all relevant time lags are considered.}
  \item{cumul}{Logical. If \code{FALSE} (the default), instantaneous path coefficients are returned, otherwise cumulative ones are provided.}
  \item{conf}{The confidence level for edges in the object of class \code{dlsem}. Default is 0.95.}
  \item{nitt}{An integer number greater or equal than 1000, representing the number of Monte Carlo replications to approximate condfidence intervals. Default is 50000.}
}
\value{A list of matrices, one for each path between starting variables and the ending variable, plus one.
Each matrix contains point estimate and \code{conf}*100\% confidence interval for the associated coefficient at each requested time lag.
The last matrix contains point estimate and \code{conf}*100\% confidence interval for the overall coefficient at each requested time lag.}
%\references{}
\author{Alessandro Magrini <magrini@disia.unifi.it>}
\seealso{\link{dlsem}}
\examples{
data(agres)
mycode <- list(
  GVA~quec(NPATENT,1,15),
  PPI~quec(NPATENT,0,9)+quec(GVA,0,10),
  ENTR_INCOME~quec(NPATENT,1,13)+quec(GVA,1,14)
  )
myfit <- dlsem(mycode,group="COUNTRY",context=c("GDP","FARM_SIZE"),
  data=agres,uniroot.check=TRUE,imputation=FALSE,log=TRUE)
pathAnal(myfit,from="NPATENT",to="PPI",lag=c(0,5,10,15),cumul=TRUE,nitt=10000)
}
