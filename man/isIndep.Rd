\name{isIndep}
\alias{isIndep}
\title{Check conditional independence}
\description{Check conditional independence between two variables using the d-separation criterion (Pearl, 2000).}
\usage{isIndep(x, var1, var2, given = NULL, conf = 0.95)}
\arguments{
  \item{x}{An object of class \code{dlsem}.}
  \item{var1}{The name of the first variable.}
  \item{var2}{The name of the second variable.}
  \item{given}{A vector containing the names of conditioning variables. If \code{NULL}, marginal independence is checked.}
  \item{conf}{The confidence level for edges in the object of class \code{dlsem}. Default is 0.95.}
}
\value{Logical}
\note{The result is unchanged if arguments \code{var1} and \code{var2} are switched.}
\references{J. Pearl (2000). Causality: Models, Reasoning, and Inference. Cambridge University Press. Cambridge, UK.}
%\author{Alessandro Magrini <magrini@disia.unifi.it>}
\seealso{\link{dlsem}}
\examples{
data(agres)
mycode <- list(
  GVA~quec(NPATENT,1,15),
  PPI~quec(NPATENT,0,13)+quec(GVA,0,14),
  ENTR_INCOME~quec(NPATENT,0,14)+quec(GVA,1,14)
  )
myfit <- dlsem(mycode,group="COUNTRY",exogenous=c("GDP","FARM_SIZE"),
  data=agres,uniroot.check=TRUE,imputation=FALSE,log=TRUE)
isIndep(myfit,"NPATENT","ENTR_INCOME",given=c("GVA","PPI"))
}
