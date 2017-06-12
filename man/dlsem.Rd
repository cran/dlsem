\name{dlsem}
\alias{dlsem}
\title{Distributed-lag linear structural equation modelling}
\description{Estimation of a distributed-lag linear structural equation model.}
\usage{dlsem(model.code, group = NULL, exogenous = NULL, data, log = FALSE, control = NULL,
  test = "adf", combine = "choi", k = 0, lshort = TRUE, maxdiff = 3, tol = 0.0001,
  maxiter = 500, selection = "aic")}
\arguments{
  \item{model.code}{A list of objects of class \code{formula}, each describing a single regression model. See \code{Details}.}
  \item{group}{The name of the group factor (optional). If \code{NULL}, no groups are considered.}
  \item{exogenous}{The name of exogenous variables (optional). Exogenous variables can be either quantitative or qualitative, must appear in no equation, and are not lagged.}  
  \item{data}{An object of class \code{data.frame} containing data.}
  \item{log}{Logical. If \code{TRUE}, logarithmic transformation is applied to strictly positive quantitative variables. Default is \code{FALSE}.}
  \item{control}{A list containing options for estimation. See \code{Details}.}
  \item{test}{The unit root test to use, that can be either \code{"adf"} or \code{"kpss"} (see \link{unirootTest}). Default is \code{"adf"}.}
  \item{combine}{The method to combine p-values of different groups, that can be either \code{"choi"} or \code{"demetrescu"} (see \link{unirootTest}).
    Ignored if \code{group} is \code{NULL}. Default is \code{"choi"}.}
  \item{k}{The lag order to calculate the statistic of the Augmented Dickey-Fuller test.
    Ignored if \code{test}=\code{"kpss"}. Default is 0.}
  \item{lshort}{Logical. If \code{TRUE}, the short version of the truncation lag parameter is used for the KPSS test.
    Ignored if \code{test}=\code{"adf"}. Default is \code{TRUE}.}
  \item{maxdiff}{The maximum differentiation order to apply. If \code{maxdiff}=0, differentiation will not be applied. Default is 3.}
  \item{maxiter}{The maximum number of iterations for the EM algorithm. Default is 500. If \code{maxiter}=0, imputation will not be performed.}
  \item{tol}{The tolerance threshold of the EM algorithm. Default is 0.0001.}
  \item{selection}{The criterion to be used for the adaptation of lag shapes, that can be one among \code{"aic"} to minimise the Akaike Information Criterion (Akaike, 1974),
    \code{"bic"} to minimise the Bayesian Information Criterion (Schwarz, 1978), and \code{"mdl"} to minimise the Minimum Description Length (Rissanen, 1978). Default is \code{"aic"}.}
}
\details{Formulas cannot contain qualitative variables or interaction terms (no ':' or '*' symbols), and
may contain the following operators for lag specification:
  \itemize{
  \item{\code{quec}: }{quadratic (2nd order polynomial) lag shape with endpoint constraints;}
  \item{\code{qdec}: }{quadratic (2nd order polynomial) decreasing lag shape;}
  \item{\code{gamma}: }{gamma lag shape.}
  }
Each operator must have the following three arguments (provided within brackets):
  \enumerate{
  \item{the name of the covariate to which the lag is applied;}
  \item{the minimum lag with a non-zero coefficient (for 2nd order polynomial lag shapes), or the \code{delta} parameter (for the gamma lag shape);}
  \item{the maximum lag with a non-zero coefficient (for 2nd order polynomial lag shapes), or the \code{lambda} parameter (for the gamma lag shape).}
  }
For example, \code{quec(X1,3,15)} indicates that a quadratic lag shape with endpoint constraints must be applied to variable X1 in the interval (3,15),
and \code{gamma(X1,0.75,0.8)} indicates that a gamma lag shape with \code{delta}=0.75 and \code{lambda}=0.8 must be applied to variable X1.
See Judge et al. (1985, Chapters 9-10) for more details.

The formula of regression models with no covariates excepting exogenous variables can be omitted from argument \code{model.code}.
The group factor and exogenous variables must not appear in any formula.

Argument \code{control} must be a named list containing one or more among the following components:
  \itemize{
  %\item{\code{L}: }{a named vector of non-negative integer values including the highest lag with non-zero autocorrelation for one or more response variables.
  %If greater than 0, the Newey-West correction of the covariance matrix of estimates (Newey and West, 1987) is used. Default is 0 for all response variables.}
  \item{\code{adapt}: }{a named vector of logical values indicating if adaptation of lag shapes must be performed for one or more response variables. Default is \code{FALSE} for all response variables.}
  \item{\code{max.gestation}: }{a named list. Each component of the list must refer to one response variable and contain a named vector, including the maximum gestation lag for one or more covariates.
  If not provided, it is taken as equal to \code{max.lead} (see below). Ignored if \code{adapt}=\code{FALSE} for a certain covariate.}
  \item{\code{max.lead}: }{a named list. Each component of the list must refer to one response variable and contain a named vector, including the maximum lead lag for one or more covariates.
  If not provided, it is computed accordingly to the sample size. Ignored if \code{adapt}=\code{FALSE} for a certain covariate.}
  \item{\code{min.width}: }{a named list. Each component of the list must refer to one response variable and contain a named vector, including the minimum lag width for one or more covariates.
  It cannot be greater than \code{max.lead}. If not provided, it is taken as 0. Ignored if \code{adapt}=\code{FALSE} for a certain covariate.}
  \item{\code{sign}: }{a named list. Each component of the list must refer to one response variable and contain a named vector, including the sign
  (either '+' for non-negative, or '-' for non-positive) of the coefficients of one or more covariates.
  If not provided, adaptation will disregard the sign of coefficients. Ignored if \code{adapt}=\code{FALSE} for a certain covariate.}
  }
}
\value{An object of class \code{dlsem}, with the following components:
  \item{estimate}{A list of objects of class \code{lm}, one for each response variable.}
  \item{model.code}{The model code after eventual adaptation.}
  \item{exogenous}{The names of exogenous variables.}
  \item{group}{The name of the group factor. \code{NULL} is returned if \code{group}=\code{NULL}.}
  \item{log}{The value provided to argument \code{log}.}
  \item{ndiff}{The order of differentiation.}
  \item{data.orig}{The dataset provided to argument \code{data}.}
  \item{data.used}{Data used in the estimation, that is after eventual logarithmic transformation and differentiation.}
S3 methods available for class \code{dlsem} are:
  \item{print}{provides essential information on the structural model.}
  \item{summary}{shows summaries of estimation.}
  \item{plot}{displays the directed acyclic graph. 
  If option \code{show.sign} is set to \code{TRUE} (the default), each significant edge is coloured with respect to the sign of its causal effect (green: positive, red: negative).
  If option \code{show.ns} is set to \code{TRUE} (the default), not statistically significant edges are shown in grey, otherwise they are omitted.}
  \item{residuals}{returns residuals.}
}
\references{
%W. K. Newey, and K. D. West (1978). A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix. \emph{Econometrica}, 55(3), 703-708.}
%\author{Alessandro Magrini <alessandro.magrini@for.unipi.it>}

H. Akaike (1974). A New Look at the Statistical Identification Model. \emph{IEEE Transactions on Automatic Control}, 19, 716-723.

G. G. Judge, W. E. Griffiths, R. C. Hill, H. Lutkepohl, and T. C. Lee (1985). The Theory and Practice of Econometrics. John Wiley & Sons, 2nd ed., New York, US-NY.

J. Rissanen (1978). Modeling by Shortest Data Description. \emph{Automatica}, 14(5): 465-658.

G. Schwarz (1978). Estimating the Dimension of a Model. \emph{Annals of Statistics}, 6, 461-464.

%P. Schmidt (1974). An Argument for the Usefulness of the Gamma Distributed Lag Model. \emph{International Economic Review}, 15(1).
}
\seealso{\link{unirootTest}; \link{causalEff}.}
\examples{
data(industry)

# estimation without control options
mycode <- list(
  Consum~quec(Job,0,5),
  Pollution~quec(Job,1,8)+quec(Consum,1,6)
  )
myfit <- dlsem(mycode,group="Region",exogenous=c("Population","GDP"),
  data=industry,log=TRUE)


### adaptation of lag shapes (takes some seconds more)
#
#mycode <- list(
#  Consum~quec(Job,0,15),
#  Pollution~quec(Job,0,15)+quec(Consum,0,15)
#  )
#                      
#mycontrol <- list(
#  adapt=c(Consum=T,Pollution=T),
#  max.gestation=list(Consum=c(Job=3),Pollution=c(Consum=3,Job=3)),
#  max.lead=list(Consum=c(Job=15),Pollution=c(Consum=15,Job=15)),
#  min.width=list(Consum=c(Job=5),Pollution=c(Consum=5,Job=5)),
#  sign=list(Consum=c(Job="+"),Pollution=c(Consum="+",Job="+"))
#  )
#
#myfit <- dlsem(mycode,group="Region",exogenous=c("Population","GDP"),data=industry,
#  control=mycontrol,log=TRUE)


# add a qualitative exogenous variable
industry[,"Policy"] <- factor(1*(industry[,"Year"]>=2006))
levels(industry[,"Policy"]) <- c("no","yes")
myfit <- dlsem(mycode,group="Region",exogenous=c("Population","GDP","Policy"),
  data=industry,log=TRUE)
  
# summaries of estimation
summary(myfit)

# directed acyclic graph
plot(myfit)

# directed acyclic graph including only statistically significant edges
plot(myfit,show.ns=FALSE)
}
