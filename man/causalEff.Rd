\name{causalEff}
\alias{causalEff}
\title{Assessment of causal effects at different time lags}
\description{Given one or more starting variables and one ending variable,
return the coefficient associated to each directed path connecting them
and the overall one at one or more time lags.}
\usage{causalEff(x, from = NULL, to = NULL, lag = NULL, cumul = FALSE, conf = 0.95,
  use.ns = FALSE)}
\arguments{
  \item{x}{An object of class \code{dlsem}.}
  \item{from}{The name of the starting variable, or a vector of names of starting variables.}
  \item{to}{The name of the ending variable.}
  \item{lag}{A non-negative integer or a vector of non-negative integers indicating the time lags at which coefficients must be computed. If \code{NULL}, all relevant time lags are considered.}
  \item{cumul}{Logical. If \code{FALSE} (the default), instantaneous coefficients are returned, otherwise cumulative ones are provided.}
  \item{conf}{The confidence level. Default is 0.95.}
  \item{use.ns}{A logical value indicating whether edges without statistically significant causal effect (at level \code{conf}) should be considered or not. If \code{FALSE} (the default), they will be ignored.}
}
\value{A list containing one matrix for each directed path between the starting variables and the ending variable and one matrix for the overall coefficient. %, and a character vector.
Each matrix contains point estimate and confidence interval (at level \code{conf}) for the associated coefficient at each requested time lag.}
%The character vector contains the names of the variables counfounding the relationship, thus the overall coefficient at each lag should be interpreted at constant values of those variables.
\note{This is an extension of causal inference described in Pearl (2012) to take into account time lags.
%trace rules for classical path analysis (Wright, 1934; Pearl, 2012).
}
\references{
J. Pearl (2012). The Causal Foundations of Structural Equation Modelling. In: R. H. Hoyle (ed.), Handbook of Structural Equation Modelling, Chapter 5. Guilford Press, New York, US-NY. ISBN: 978-16-062-3077-0
%
%S. Wright (1934). The Method of Path Coefficients. \emph{Annals of Mathematical Statistics}, 5(3): 161-215.
}
\seealso{\link{dlsem}; \link{lagPlot}.}
\examples{
data(industry)
mycode <- list(
  Consum~quec(Job,0,5),
  Pollution~quec(Job,1,8)+quec(Consum,1,6)
  )
myfit <- dlsem(mycode,group="Region",exogenous=c("Population","GDP"),data=industry,log=TRUE)
causalEff(myfit,from="Job",to="Pollution",lag=c(0,5,10,15),cumul=TRUE)
}