\name{lagShapes}
\alias{lagShapes}
\title{Estimated lag shapes}
\description{Estimated lag shapes and their standard errors are provided.}
\usage{lagShapes(x, cumul = FALSE)}
\arguments{
  \item{x}{An object of class \code{dlsem}.}
  \item{cumul}{Logical. If \code{TRUE}, cumulative causal effects are returned. Default is \code{FALSE}.}
}
\value{A list of lists, one for each endogenous variable, each containing several matrices including the estimated lag shapes and their standard errors.
}
%\note{}
%\references{}
\seealso{\link{dlsem}; \link{causalEff}; \link{lagPlot}.}
\examples{
data(industry)
indus.code <- list(
  Consum~quec.lag(Job,0,5),
  Pollution~quec.lag(Job,2,9)+quec.lag(Consum,1,7)
  )
indus.mod <- dlsem(indus.code,group="Region",exogenous=c("Population","GDP"),data=industry,
  log=TRUE)
lagShapes(indus.mod)
}
