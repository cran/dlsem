\name{quec.lag}
\alias{quec.lag}
\alias{qdec.lag}
\alias{gamm.lag}
\title{Lag shape constructors}
\description{Lag shape constructors to be used in model formulas.}
\usage{quec.lag(x, a, b, x.group = NULL, nlag=NULL)
qdec.lag(x, a, b, x.group = NULL, nlag=NULL)
gamm.lag(x, delta, lambda, x.group = NULL, nlag=NULL)
}
\arguments{
  \item{x}{The name of the variable.}
  \item{a}{The gestation lag.}
  \item{b}{The lead lag.}
  \item{delta, lambda}{Delta and lambda parameters for the gamma lag shape (only for \code{gamm.lag}).}
  \item{nlag}{The number of lags considered. If \code{NULL} or less than the lead lag, it is set equal to the lead lag.}
  \item{x.group}{The name of the group factor (optional).}
}
%\value{}
\note{The constructor \code{gamma.lag()} is deprecated and it will be replaced by \code{gamm.lag()}.
%The lead lag of a gamma lag shape is determined numerically.
}
%\references{}
%\seealso{\link{dlsem}.}
\examples{
data(industry)
# example in linear regression
m0 <- lm(Consum ~ -1+Region+quec.lag(Job,0,5,x.group=Region), data=industry)
m1 <- lm(Consum ~ -1+Region+qdec.lag(Job,0,5,x.group=Region), data=industry)
m2 <- lm(Consum ~ -1+Region+gamm.lag(Job,0.5,0.5,x.group=Region), data=industry)
}
