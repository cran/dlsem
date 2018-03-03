\name{auto.lagPlot}
\alias{auto.lagPlot}
\title{Automated plot of lag shapes}
\description{All the single-edge pathwise causal lag shapes are saved as pdf files.}
\usage{auto.lagPlot(x, cumul = FALSE, conf = 0.95, plotDir = NULL)}
\arguments{
  \item{x}{An object of class \code{dlsem}.}
  \item{cumul}{Logical. If \code{FALSE} (the default), instantaneous causal effects are displayed, otherwise cumulative ones are provided.}
  \item{conf}{The confidence level for each plot. Default is 0.95.}
  \item{plotDir}{The directory where to save the plots. If \code{NULL} (the default), plots will be saved in the current working directory.}
}
%\value{}
%\references{}
\seealso{\link{dlsem}; \link{lagPlot}.}
\examples{
data(industry)
indus.code <- list(
  Consum~quec.lag(Job,0,4),
  Pollution~quec.lag(Job,2,6)+quec.lag(Consum,1,5)
  )
indus.mod <- dlsem(indus.code,group="Region",exogenous=c("Population","GDP"),data=industry,
  log=TRUE)
auto.lagPlot(indus.mod,plotDir=getwd())
}
