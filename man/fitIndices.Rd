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
The second list contains the Akaike Information Criterion for each regression model and the overall one.
The third list contains the Bayesian Information Criterion for each regression model and the overall one.
The fourth list contains the Minimum Information Length for each regression model and the overall one.
}
%\note{}
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
fitIndices(myfit)
}
