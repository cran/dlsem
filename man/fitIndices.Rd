\name{fitIndices}
\alias{fitIndices}
\title{Fit indices}
\description{Compute local and global fit indices.}
\usage{fitIndices(x)}
\arguments{
  \item{x}{An object of class \code{dlsem}.}
  }
\value{A list of three numeric vectors.
The fist list contains the R squared for each regression model and the overall one.
The second list contains the Akaike Information Criterion (AIC) for each regression model and the overall one.
The third list contains the Bayesian Information Criterion (BIC) for each regression model and the overall one.
}
%\note{}
%\references{}
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
fitIndices(myfit)
}
