\name{compareModels}
\alias{compareModels}
\title{Comparison among different distributed-lag linear structural equation models}
\description{Several competing distributed-lag linear structural equation models are compared based on information criteria.}
\usage{compareModels(x)}
\arguments{
  \item{x}{A list of 2 or more objects of class \code{dlsem} estimated on the same data.}
}
\value{A data.frame with one record for each model in \code{x} on the following quantities: log-likelihood, number of parameters, Akaike Information Criterion (AIC), corrected Akaike Information Criterion (AICc) and Bayesian Information criterion (BIC).}
\note{In order to keep the data constant, all the models are estimated from scratch after applying a fixed number of lags to the data, equal to the hightest lead lag across all the models.}
\seealso{\link{dlsem}.}
\examples{
data(industry)

# model with endpoint-contrained quadratic lag shapes
indus.code <- list(
  Consum~quec.lag(Job,0,5),
  Pollution~quec.lag(Job,2,9)+quec.lag(Consum,1,7)
  )
indus.mod <- dlsem(indus.code,group="Region",exogenous=c("Population","GDP"),data=industry,
  log=TRUE)
  
# model with gamma lag shapes
indus.code_2 <- list(
  Consum~gamm.lag(Job,0.7,0.05),
  Pollution~gamm.lag(Job,0.95,0.05)+gamm.lag(Consum,0.9,0.15)
  )
indus.mod_2 <- dlsem(indus.code_2,group="Region",exogenous=c("Population","GDP"),data=industry,
  log=TRUE)
  
compareModels(list(indus.mod,indus.mod_2))
}