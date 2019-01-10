\name{dlsem-package}
\alias{dlsem-package}
\docType{package}
\title{Distributed-lag linear structural equation models}
\description{
Inference functionalities for distributed-lag linear structural equation models (DLSEMs).
DLSEMs are Markovian structural causal models where each factor of the joint probability distribution
is a distributed-lag linear regression model (Magrini, 2018).
DLSEMs account for temporal delays in the dependence relationships among the variables and
allow to perform dynamic causal inference by assessing causal effects at different time lags.
Endpoint-constrained quadratic, quadratic decreasing and gamma lag shapes are available.
The main functions of the package are:
\itemize{
  \item{\link{dlsem}, }{to perform parameter estimation;}
  \item{\link{causalEff}, }{to compute all the pathwise causal lag shapes and the overall one connecting two or more variables;}
  \item{\link{lagPlot}, }{to display a pathwise or an overall causal lag shape.}
  }
}
\details{
  \tabular{ll}{
  Package: \tab dlsem\cr
  Type: \tab Package\cr
  Version: \tab 2.4.1\cr
  Date: \tab 2019-01-11\cr
  License: \tab GPL-2
  }
}
\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
A. Magrini (2018). Linear Markovian Models for Lag Exposure Assessment. Biometrical Letters, 55(2): 179-195. DOI: 10.2478/bile-2018-0012.

%G. G. Judge, W. E. Griffiths, R. C. Hill, H. Lutkepohl, and T. C. Lee (1985). The Theory and Practice of Econometrics. John Wiley \& Sons, 2nd ed., New York, US-NY. ISBN: 978-0-471-89530-5

%J. Pearl (2000). Causality: Models, Reasoning, and Inference. Cambridge University Press. Cambridge, UK. ISBN: 978-0-521-89560-6
}        
