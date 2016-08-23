\name{pathAnal}
\alias{pathAnal}
\title{Perform path analysis at different time lags}
\description{Compute the coefficient associated to each path between one or more starting variables and one ending variable at one or more time lags.}
\usage{pathAnal(x, from = NULL, to = NULL, lag = NULL, cumul = FALSE, conf = 0.95)}
\arguments{
  \item{x}{An object of class \code{dlsem}.}
  \item{from}{The name of the starting variable, or a vector of names of starting variables.}
  \item{to}{The name of the ending variable.}
  \item{lag}{A non-negative integer or a vector of non-negative integers indicating the time lag at which path coefficients must be computed. If \code{NULL}, all relevant time lags are considered.}
  \item{cumul}{Logical. If \code{FALSE} (the default), instantaneous path coefficients are returned, otherwise cumulative ones are provided.}
  \item{conf}{The confidence level for edges in the object of class \code{dlsem}. Default is 0.95.}
}
\value{A list containing one matrix for each path between the starting variables and the ending variable and one matrix for the overall coefficient. %, and a character vector.
Each matrix contains point estimate and asymptotic \code{conf}*100\% confidence interval for the associated coefficient at each requested time lag.}
%The character vector contains the names of the variables counfounding the relationship, thus the overall coefficient at each lag should be interpreted at constant values of those variables.
\note{This is an extension of trace rules for classical path analysis (Wright, 1934; Pearl, 2012).
%to distributed-lag structural equation models with second-order polynomial lag shapes.
}
\references{
J. Pearl (2012). The Causal Foundations of Structural Equation Modelling. In: R. H. Hoyle (ed.), Handbook of Structural Equation Modelling, Chapter 5. Guilford Press, New York, US-NY.

%J. Pearl (2009). Causal Inference in Statistics: An overview. \emph{Statistics Surveys}, 3: 96-146.

S. Wright (1934). The Method of Path Coefficients. \emph{Annals of Mathematical Statistics}, 5(3): 161-215.
}
\author{Alessandro Magrini <magrini@disia.unifi.it>}
\seealso{\link{dlsem}}
\examples{
data(agres)
mycode <- list(
  GVA~quec(NPATENT,1,15),
  PPI~quec(NPATENT,0,13)+quec(GVA,0,14),
  ENTR_INCOME~quec(NPATENT,0,14)+quec(GVA,1,14)
  )
myfit <- dlsem(mycode,group="COUNTRY",context=c("GDP","FARM_SIZE"),
  data=agres,uniroot.check=TRUE,imputation=FALSE,log=TRUE)
pathAnal(myfit,from="NPATENT",to="PPI",lag=c(0,5,10,15),cumul=TRUE)
}
