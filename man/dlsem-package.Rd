\name{dlsem-package}
\alias{dlsem-package}
\docType{package}
\title{Distributed-Lag Structural Equation Modelling}
\description{Fit distributed-lag structural equation models and perform path analysis at different time lags.
%A distributed-lag structural equation model can be fitted using function \link{dlsem}.
%Path analysis can be performed using function \link{pathAnal}.
%Function \link{lagPlot} displays the lag shape associated to a causal effect or a path.
}
\details{
  \tabular{ll}{
  Package: \tab dlsem\cr
  Type: \tab Package\cr
  Version: \tab 1.0\cr
  Date: \tab 2016-05-19\cr
  License: \tab GPL-2
  }
}
\author{Alessandro Magrini <magrini@disia.unifi.it>}
\references{A. Magrini, F. Bartolini, A. Coli, and B. Pacini (2016). Distributed-Lag Structural Equation Modelling:
An Application to Impact Assessment of Research Activity on European Agriculture.
\emph{Proceedings of the 48th Meeting of the Italian Statistical Society}, 8-10 June 2016, Salerno, IT.}
%\examples{
%# load data on agricultural research
%data(agres)
%
%# model code
%mycode <- list(
%  GVA~quec(NPATENT,0,10),
%  PPI~quec(NPATENT,0,10)+quec(GVA,0,10),
%  ENTR_INCOME~quec(NPATENT,0,10)+quec(GVA,0,10)
%  )
%
%# define control options
%mycontrol <- list(
%  adapt=c(GVA=TRUE,PPI=TRUE,ENTR_INCOME=TRUE),
%  max.gestation=list(GVA=c(NPATENT=3),PPI=c(NPATENT=3,GVA=3),ENTR_INCOME=c(NPATENT=3,GVA=3)),
%  min.width=list(GVA=c(NPATENT=3),PPI=c(NPATENT=3,GVA=3),ENTR_INCOME=c(NPATENT=3,GVA=3)),
%  sign=list(GVA=c(NPATENT="+"),PPI=c(NPATENT="-",GVA="-"),ENTR_INCOME=c(NPATENT="+",GVA="+"))
%  )
%
%# fit the dlsem
%mod0 <- dlsem(mycode,group="COUNTRY",context=c("GDP","FARM_SIZE"),data=agres,
%  control=mycontrol,uniroot.check=TRUE,maxdiff=1,imputation=FALSE,log=TRUE)
%
%# see summaries
%summary(mod0)
%
%# plot the directed acyclic graph of the model with significant edges
%plot(mod0)
%
%# path analysis example
%#   at first 10 time lags
%pathAnal(mod0,from="NPATENT",to="ENTR_INCOME",lag=0:10,nitt=10000)
%#   cumulative after 10 time lags
%pathAnal(mod0,from="NPATENT",to="ENTR_INCOME",lag=10,cumul=TRUE,nitt=10000)
%
%# plot the lag shape associated to specific paths
%lagPlot(mod0,path="NPATENT*PPI",nitt=10000)
%lagPlot(mod0,path="NPATENT*GVA*PPI",nitt=10000)
%
%# plot the lag shape associated to an overall causal effect
%lagPlot(mod0,from="NPATENT",to="PPI",nitt=10000)
%}        
