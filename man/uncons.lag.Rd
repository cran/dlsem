\name{uncons.lag}
\alias{quec.lag}
\alias{qdec.lag}
\alias{gamm.lag}
\alias{uncons.lag}
\alias{almon.lag}
\title{Lag shape transformation for model formulas}
\description{Lag shape transformation to be used in model formulas.}
\usage{quec.lag(x, a, b, x.group = NULL)
qdec.lag(x, a, b, x.group = NULL)
gamm.lag(x, delta, lambda, x.group = NULL)
uncons.lag(x, a, b, x.group = NULL)
almon.lag(x, a, b, p, x.group = NULL)
}
\arguments{
  \item{x}{The name of the variable.}
  \item{a}{The gestation lag.}
  \item{b}{The lead lag.}
  \item{p}{The polynomial order for the Almon's lag shape (only for \code{almon.lag}).}
  \item{delta, lambda}{Delta and lambda parameters for the gamma lag shape (only for \code{gamm.lag}).}
  \item{x.group}{The name of the group factor (optional).}
}
%\value{}
%\note{}
%\references{}
%\seealso{\link{dlsem}.}
\examples{
data(industry)
# example in linear regression
m0 <- lm(Consum ~ -1+Region+quec.lag(Job,0,5,x.group=Region), data=industry)
m1 <- lm(Consum ~ -1+Region+qdec.lag(Job,0,5,x.group=Region), data=industry)
m2 <- lm(Consum ~ -1+Region+gamm.lag(Job,0.5,0.5,x.group=Region), data=industry)
m3 <- lm(Consum ~ -1+Region+uncons.lag(Job,0,3,x.group=Region), data=industry)
m4 <- lm(Consum ~ -1+Region+almon.lag(Job,0,3,2,x.group=Region), data=industry)
}
