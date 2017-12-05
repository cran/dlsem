\name{dlsem-package}
\alias{dlsem-package}
\docType{package}
\title{Distributed-lag linear structural equation models}
\description{
Inference functionalities for distributed-lag linear structural equation models (DLSEMs).
DLSEMs were recently proposed by Magrini et al. (2016) as an extension of
linear structural causal models (Pearl, 2000, Chpater 5), where each factor of the
joint probability distribution is a distributed-lag linear regression with
constrained lag shapes (Judge et al., 1985, Chapters 9-10).
DLSEMs account for temporal delays in the dependence relationships among the variables
and allow to assess causal effects at different time lags.
Endpoint-constrained quadratic, quadratic decreasing and gamma lag shapes are available.
The main functions of the package are:
\itemize{
  \item{\link{dlsem}, }{to perform parameter estimation;}
  \item{\link{causalEff}, }{to assess causal effects at different time lags;}
  \item{\link{lagPlot}, }{to plot the lag shape associated to a path or to an overall causal effect.} 
  }
}
\details{
  \tabular{ll}{
  Package: \tab dlsem\cr
  Type: \tab Package\cr
  Version: \tab 2.0\cr
  Date: \tab 2017-12-05\cr
  License: \tab GPL-2
  }
}
\author{Alessandro Magrini <magrini@disia.unifi.it>}
\references{
G. G. Judge, W. E. Griffiths, R. C. Hill, H. Lutkepohl, and T. C. Lee (1985). The Theory and Practice of Econometrics. John Wiley \& Sons, 2nd ed., New York, US-NY. ISBN: 978-0-471-89530-5

A. Magrini, F. Bartolini, A. Coli, and B. Pacini (2016). Distributed-Lag Structural Equation Modelling:
An Application to Impact Assessment of Research Activity on European Agriculture.
\emph{Proceedings of the 48th Meeting of the Italian Statistical Society}, 8-10 June 2016, Salerno, IT. ISBN: 978-88-619-7061-8

J. Pearl (2000). Causality: Models, Reasoning, and Inference. Cambridge University Press. Cambridge, UK. ISBN: 978-0-521-89560-6
}        
