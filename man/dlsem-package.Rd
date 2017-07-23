\name{dlsem-package}
\alias{dlsem-package}
\docType{package}
\title{Distributed-lag linear structural equation modelling}
\description{Inference functionalities for distributed-lag linear structural equation models (Magrini et al., 2016).
A distributed-lag linear structural equation model is a structural equation model (Pearl, 2000, Chapter 5)
composed of distributed-lag linear regression models (Judge, 1985, Chapters 9-10).
Endpoint-constrained quadratic, quadratic decreasing and gamma lag shapes are available.
The main functions of the package are:
\itemize{
  \item{\link{dlsem}, }{to perform parameter estimation;}
  \item{\link{causalEff}, }{to assess %and decompose
  causal effects at different time lags;}
  \item{\link{lagPlot}, }{to plot the lag shape associated to a path or to an overall causal effect.} 
  }
}
\details{
  \tabular{ll}{
  Package: \tab dlsem\cr
  Type: \tab Package\cr
  Version: \tab 1.9\cr
  Date: \tab 2017-07-24\cr
  License: \tab GPL-2
  }
}
\author{Alessandro Magrini <alessandro.magrini@for.unipi.it>}
\references{
G. G. Judge, W. E. Griffiths, R. C. Hill, H. Lutkepohl, and T. C. Lee (1985). The Theory and Practice of Econometrics.
John Wiley \& Sons, 2nd ed., New York, US-NY.

A. Magrini, F. Bartolini, A. Coli, and B. Pacini (2016). Distributed-Lag Structural Equation Modelling:
An Application to Impact Assessment of Research Activity on European Agriculture.
\emph{Proceedings of the 48th Meeting of the Italian Statistical Society}, 8-10 June 2016, Salerno, IT.

%J. Pearl (2012). The Causal Foundations of Structural Equation Modelling.
%In: R. H. Hoyle (ed.), Handbook of Structural Equation Modelling, Chapter 5. Guilford Press, New York, US-NY. 

J. Pearl (2000). Causality: Models, Reasoning, and Inference.
Cambridge University Press. Cambridge, UK.
}        
