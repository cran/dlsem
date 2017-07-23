\name{makePrediction}
\alias{makePrediction}
\title{Make prediction}
\description{Make prediction from a distributed-lag linear structural equation model.}
\usage{makePrediction(x, newdata = NULL, group.levels = NULL, n = NULL, conf = 0.95)}
\arguments{
  \item{x}{An object of class \code{dlsem}.}
  \item{newdata}{An object of \code{data.frame} containing data on exogenous variables to condition on (optional).
  If omitted, sampling is conditioned to data on exogenous variables used in the estimation.
  Ignored if the model includes no exogenous variables.}
  \item{group.levels}{The name of groups (levels of the group factor) for which prediction should be made (optional).
  If omitted, prediction is made for the same groups as in data used in the estimation.
  Ignored if the model does not include a group factor, or if the model includes one or more exogenous variables.}
  \item{n}{The number of periods from time 0 (the first period in data used in the estimation) for which prediction should be made (optional).
  If omitted, prediction is made for the same periods as in data used in the estimation.
  Ignored if the model includes one or more exogenous variables.}
  \item{conf}{The confidence level. Default is 0.95.}
}
%\details{}
\value{A list of three objects of class \code{data.frame}.
The first contains pointwise prediction, while the second and the third contain the lower and the upper bounds of the predictive confidence interval (at level \code{conf}), respectively.
If logarithmic transformation was applied before the estimation, predicted values are on the logarithmic scale.
If k-th order differentiation was applied before the estimation, predicted values are differences with the k-th previous period.
If both logarithmic transformation and k-th order differentiation were applied before the estimation, predicted values are percentage variations with respect to the k-th previous period.
}
%\note{}
%\references{}
%\author{Alessandro Magrini <alessandro.magrini@for.unipi.it>}
\seealso{\link{dlsem}.}
\examples{
data(industry)
mycode <- list(
  Consum~quec(Job,0,5),
  Pollution~quec(Job,1,8)+quec(Consum,1,6)
  )
myfit <- dlsem(mycode,group="Region",exogenous=c("Population","GDP"),data=industry,log=TRUE)

### predicted values are percentage variations with respect to the previous period 
#
## prediction for the same data used in the estimation
makePrediction(myfit)
#
## prediction for new data on exogenous variables
newD <- data.frame(
  Region=factor(rep(5,20)),Population=rnorm(20,1e+6,1e+5),GDP=rnorm(20,1e+5,1e+4)
  )
makePrediction(myfit,newdata=newD)
}
