\name{lagPlot}
\alias{lagPlot}
\title{Plot lag shapes}
\description{Plot the lag shape associated to a path or to an overall causal effect.}
\usage{lagPlot(model, from = NULL, to = NULL, path = NULL, maxlag = NULL, cumul = FALSE,
  conf = 0.95, ylim = NULL, title = NULL)}
\arguments{
  \item{model}{An object of class \code{dlsem}.}
  \item{from}{The name of the covariate.}
  \item{to}{The name of the response variable.}
  \item{path}{The name of the path, indicated as a string made of the names of the variables in the path separated by '*'.
  Ignored if both \code{from} and \code{to} are not \code{NULL}.}
  \item{maxlag}{The maximum lag displayed (optional).}
  \item{cumul}{Logical. If \code{FALSE} (the default), instantaneous coefficients are returned, otherwise cumulative ones are provided.}
  \item{conf}{The confidence level for edges in the object of class \code{dlsem}. Default is 0.95.}
  \item{ylim}{A vector of two numerical values indicating the limits of the y axis (optional). If \code{NULL}, the limits of the y axis are computed automatically.}
  \item{title}{The title of the plot (optional). If \code{NULL}, a default title is used.}
}
%\value{}
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

# the lag shape associated to specific paths
lagPlot(myfit,path="Job*Pollution")
lagPlot(myfit,path="Job*Consum*Pollution")

# the lag shape associated to an overall causal effect
lagPlot(myfit,from="Job",to="Pollution")
}
