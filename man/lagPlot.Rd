\name{lagPlot}
\alias{lagPlot}
\title{Plot lag shapes}
\description{Plot the lag shape associated to a causal effect or a path.}
\usage{lagPlot(model, from = NULL, to = NULL, path = NULL, maxlag = NULL, cumul = FALSE,
  conf = 0.95, nitt = 50000, ylim = NULL, title = NULL)}
\arguments{
  \item{model}{An object of class \code{lm} or \code{dlsem}.}
  \item{from}{The name of the covariate. Ignored if an object of class \code{lm} is provided to argument \code{model}.}
  \item{to}{The name of the response variable.}
  \item{path}{The name of the path. Ignored if an object of class \code{lm} is provided to argument \code{model}, or if both \code{from} and \code{to} are not \code{NULL}.}
  \item{maxlag}{The maximum lag displayed (optional).}
  \item{cumul}{Logical. If \code{FALSE} (the default), instantaneous path coefficients are returned, otherwise cumulative ones are provided.}
  \item{conf}{The confidence level for edges in the object of class \code{dlsem}. Default is 0.95.}
  \item{nitt}{An integer number greater or equal than 1000, representing the number of Monte Carlo replications to approximate confidence intervals. Default is 50000.}
  \item{ylim}{A vector of two numerical values indicating the limits of the y axis (optional). If \code{NULL}, the limits of the y axis are computed automatically.}
  \item{title}{The title of the plot (optional). If \code{NULL}, a default title is used.}
}
%\value{}
%\references{}
\author{Alessandro Magrini <magrini@disia.unifi.it>}
\seealso{\link{dlaglm}, \link{dlsem}}
\examples{
data(agres)
mycode <- list(
  GVA~quec(NPATENT,0,4),
  PPI~quec(NPATENT,0,6)+quec(GVA,0,10),
  ENTR_INCOME~quec(NPATENT,0,3)+quec(GVA,1,10)
  )
mod0 <- dlsem(mycode,group="COUNTRY",context=c("GDP","FARM_SIZE"),
  data=agres,uniroot.check=TRUE,imputation=FALSE,log=TRUE)

# the lag shape associated to specific paths
lagPlot(mod0,path="NPATENT*PPI",nitt=10000)
lagPlot(mod0,path="NPATENT*GVA*PPI",nitt=10000)

# the lag shape associated to an overall causal effect
lagPlot(mod0,from="NPATENT",to="PPI",nitt=10000)
}
