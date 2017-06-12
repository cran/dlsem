\name{isIndep}
\alias{isIndep}
\title{Check conditional independence}
\description{Check conditional independence between two variables using the d-separation criterion (Pearl, 2000).}
\usage{isIndep(x, var1 = NULL, var2 = NULL, given = NULL, conf = 0.95, use.ns = FALSE)}
\arguments{
  \item{x}{An object of class \code{dlsem}.}
  \item{var1}{The name of the first variable.}
  \item{var2}{The name of the second variable.}
  \item{given}{A vector containing the names of conditioning variables. If \code{NULL}, marginal independence is checked.}
  \item{conf}{The confidence level for each edge: only statistically significant edges at such level are considered. Default is 0.95.}
  \item{use.ns}{A logical value indicating whether not statistically significant edges (at level \code{conf}) should be considered or not. If \code{FALSE} (the default), they will be ignored.}
}
\value{Logical}
\note{The result is unchanged if arguments \code{var1} and \code{var2} are switched.}
\references{J. Pearl (2000). Causality: Models, Reasoning, and Inference. Cambridge University Press. Cambridge, UK.}
%\author{Alessandro Magrini <alessandro.magrini@for.unipi.it>}
\seealso{\link{dlsem}.}
\examples{
data(industry)
mycode <- list(
  Consum~quec(Job,0,5),
  Pollution~quec(Job,1,8)+quec(Consum,1,6)
  )
myfit <- dlsem(mycode,group="Region",exogenous=c("Population","GDP"),
  data=industry,log=TRUE)
isIndep(myfit,"Job","Pollution",given=c("Consum"))
}
