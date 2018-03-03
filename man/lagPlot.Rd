\name{lagPlot}
\alias{lagPlot}
\title{Plot lag shapes}
\description{A pathwise or an overall causal lag shape is displayed.}
\usage{lagPlot(x, from = NULL, to = NULL, path = NULL, maxlag = NULL, cumul = FALSE,
  conf = 0.95, use.ns = FALSE, ylim = NULL, title = NULL)}
\arguments{
  \item{x}{An object of class \code{dlsem}.}
  \item{from, to, path}{
  To display the overall causal lag shape of a variable to another one,
  their names must be provided to arguments \code{from} and \code{to}, respectively.
  To display a pathwise causal lag shape,
  the name of the path, indicated as a string made of the names of the variables in the path separated by '*', must be provided to argument \code{path}.
  Argument \code{path} will be ignored if both \code{from} and \code{to} are not \code{NULL}.}
  \item{maxlag}{The maximum lag displayed (optional).}
  \item{cumul}{Logical. If \code{FALSE} (the default), instantaneous causal effects are returned, otherwise cumulative ones are provided.}
  \item{conf}{The confidence level for each edge: only statistically significant edges at such level are considered. Default is 0.95.}
  \item{use.ns}{A logical value indicating whether not statistically significant edges (at level \code{conf}) should be considered or not. If \code{FALSE} (the default), they will be ignored.}
  \item{ylim}{A vector of two numerical values indicating the limits of the y axis (optional). If \code{NULL}, the limits of the y axis are computed automatically.}
  \item{title}{The title of the plot (optional). If \code{NULL}, a default title is used.}
}
%\value{}
%\references{}
\seealso{\link{dlsem}; \link{causalEff}; \link{lagPlot}.}
\examples{
data(industry)
indus.code <- list(
  Consum~quec.lag(Job,0,4),
  Pollution~quec.lag(Job,2,6)+quec.lag(Consum,1,5)
  )
indus.mod <- dlsem(indus.code,group="Region",exogenous=c("Population","GDP"),data=industry,
  log=TRUE)

# the lag shape of the causal effect associated to specific paths
lagPlot(indus.mod,path="Job*Pollution")
lagPlot(indus.mod,path="Job*Consum*Pollution")

# the lag shape of an overall causal effect
lagPlot(indus.mod,from="Job",to="Pollution")
}
