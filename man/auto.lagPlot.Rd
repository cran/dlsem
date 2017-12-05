\name{auto.lagPlot}
\alias{auto.lagPlot}
\title{Automated plot of lag shapes}
\description{Plot the lag shape associated to each edge.}
\usage{auto.lagPlot(x, cumul = FALSE, conf = 0.95, plotDir = NULL)}
\arguments{
  \item{x}{An object of class \code{dlsem}.}
  \item{cumul}{Logical. If \code{FALSE} (the default), instantaneous coefficients are returned, otherwise cumulative ones are provided.}
  \item{conf}{The confidence level for each plot. Default is 0.95.}
  \item{plotDir}{The directory where to save the plots. If \code{NULL} (the default), plots will be saved in the current working directory.}
}
%\value{}
%\references{}
\seealso{\link{dlsem}; \link{lagPlot}.}
\examples{
data(industry)
mycode <- list(
  Consum~quec(Job,0,5),
  Pollution~quec(Job,1,8)+quec(Consum,1,6)
  )
myfit <- dlsem(mycode,group="Region",exogenous=c("Population","GDP"),data=industry,log=TRUE)
auto.lagPlot(myfit,plotDir=getwd())
}
