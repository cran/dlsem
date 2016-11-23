\name{dlsem}
\alias{dlsem}
\title{Distributed-lag structural equation modelling}
\description{Estimation of a structural equation model with second-order polynomial and gamma lag shapes.}
\usage{dlsem(model.code, group = NULL, exogenous = NULL, data, log = FALSE, control = NULL,
  imputation = TRUE, uniroot.check = TRUE, test = "adf", combine = "choi", k = 0,
  lshort = TRUE, maxdiff = 5, tol = 0.0001, maxit = 500, selection = "aic",
  plotDir = NULL)}
\arguments{
  \item{model.code}{A list of objects of class \code{formula}, each describing a single regression model. See \code{Details}.}
  \item{group}{The name of the group factor (optional). If \code{NULL}, no groups are considered.}
  \item{exogenous}{The name of exogenous variables (optional). Exogenous variables never appear on the left side of an equation and are not lagged.}  
  \item{data}{An object of class \code{data.frame} containing data.}
  \item{log}{Logical. If \code{TRUE}, logarithmic transformation is applied to quantitative variables. Default is \code{FALSE}.}
  \item{control}{A list containing options for estimation. See \code{Details}.}
  \item{imputation}{Logical. If \code{TRUE}, missing values will be imputed using the EM algorithm. Default is \code{FALSE}.}
  \item{uniroot.check}{Logical. If \code{TRUE}, unit root test is performed for each variable, and appropriate differentation is applied. Default is \code{FALSE}.}
  \item{test}{The unit root test to use, that can be either \code{"adf"} or \code{"kpss"} (see \link{unirootTest}). Ignored if \code{uniroot.check}=\code{FALSE}. Default is \code{"adf"}.}
  \item{combine}{The method to combine p-values of different groups, that can be either \code{"choi"} or \code{"demetrescu"} (see \link{unirootTest}).
    Ignored if \code{uniroot.check}=\code{FALSE} or \code{group} is \code{NULL}. Default is \code{"choi"}.}
  \item{k}{The lag order to calculate the statistic of the Augmented Dickey-Fuller test.
    Ignored if \code{uniroot.check}=\code{FALSE} or if \code{test}=\code{"kpss"}. Default is 0.}
  \item{lshort}{Logical. If \code{TRUE}, the short version of the truncation lag parameter is used for the KPSS test.
    Ignored if \code{uniroot.check}=\code{FALSE} or if \code{test}=\code{"adf"}. Default is \code{TRUE}.}
  \item{maxdiff}{The maximum differentiation order to apply. Ignored if \code{uniroot.check}=\code{FALSE}. Default is 5.}
  \item{maxit}{The maximum number of iterations for the EM algorithm (see \link{EM.imputation}). Ignored if \code{imputation}=\code{FALSE}. Default is 500.}
  \item{tol}{The tolerance threshold of the EM algorithm (see \link{EM.imputation}). Ignored if \code{imputation}=\code{FALSE}. Default is 0.0001.}
  \item{selection}{The criterion to be used for the adaptation of lag shapes, that can be one among \code{"aic"} to minimise the Akaike Information Criterion,
    \code{"bic"} to minimise the Bayesian Information Criterion (AIC), and \code{"mdl"} to minimise the Minimum Description Length. Default is \code{"aic"}.}
  %\item{prec.gamma}{A real number strictly contained in (0, 0.1) indicating the tolerance in the adaptation of gamma lag shapes.
  %  The lower the tolerance, the more the models which will be considered during adaptation of lag shapes.
  %  Since the time of completion increases exponentially in the decreasing of tolerance, tolerance values less than 0.01 are not recommended. Default is 0.025.}
  \item{plotDir}{A directory where to save the plots of the lag shapes (optional). If \code{NULL}, no plots will be produced.}
}
\details{Formulas cannot contain interaction terms (no ':' or '*' symbols), and
may contain the following operators for lag specification:
  \itemize{
  \item{\code{quec}: }{quadratic (2nd order polynomial) lag shape with endpoint constraints;}
  \item{\code{qdec}: }{quadratic (2nd order polynomial) decreasing lag shape.}
  \item{\code{gamma}: }{gamma lag shape.}
  }
Each operator must have the following three arguments (provided within brackets):
  \enumerate{
  \item{the name of the covariate to which the lag is applied;}
  \item{the minimum lag with a non-zero coefficient (for 2nd order polynomial lag shapes), or the \code{delta} parameter (for the gamma lag shape)}
  \item{the maximum lag with a non-zero coefficient (for 2nd order polynomial lag shapes), or the \code{lambda} parameter (for the gamma lag shape).}
  }
For example, \code{quec(X1,3,15)} indicates that a quadratic lag shape with endpoint constraints must be applied to variable X1 in the interval (3,15),
and \code{gamma(X1,0.75,0.8)} indicates that a gamma lag shape with \code{delta}=0.75 and \code{lambda}=0.8 must be applied to variable X1.
The formula of regression models with no covariates excepting exogenous variables can be omitted from argument \code{model.code}.
Variables appearing in any formula are treated as quantitative.
The group factor and exogenous variables must not appear in any formula.

Argument \code{control} must be a named list containing one or more among the following components:
  \itemize{
  \item{\code{L}: }{a named vector of non-negative integer values including the highest lag with non-zero autocorrelation for one or more response variables.
  If greater than 0, the Newey-West correction of the covariance matrix of estimates (Newey and West, 1987) is used. Default is 0 for all response variables.}
  \item{\code{adapt}: }{a named vector of logical values indicating if adaptation of lag shapes must be performed for one or more response variables. Default is \code{FALSE} for all response variables.}
  \item{\code{max.gestation}: }{a named list. Each component of the list must refer to one response variable and contain a named vector, including the maximum gestation lag for one or more covariates.
  If not provided, it is taken as equal to \code{max.width} (see below). Ignored if \code{adapt}=\code{FALSE} for a certain covariate.}
  \item{\code{min.width}: }{a named list. Each component of the list must refer to one response variable and contain a named vector, including the minimum lag width for one or more covariates.
  If not provided, it is taken as 0. Ignored if \code{adapt}=\code{FALSE} for a certain covariate.}
  \item{\code{max.width}: }{a named list. Each component of the list must refer to one response variable and contain a named vector, including the maximum lag width for one or more covariates.
  If not provided, it is computed accordingly to the sample size. Ignored if \code{adapt}=\code{FALSE} for a certain covariate.}
  \item{\code{sign}: }{a named list. Each component of the list must refer to one response variable and contain a named vector, including the sign
  (either '+' for non-negative, or '-' for non-positive) of the coefficients of one or more covariates.
  If not provided, adaptation will disregard the sign of coefficients. Ignored if \code{adapt}=\code{FALSE} for a certain covariate.}
  }
%For the last three components, values associated to covariates without a lag specification
%or to variables not appearing in the model code will be ignored.
Variables appearing in the model code but not included in data will be considered as unobserved.
If there is at least one unobserved variable, imputation using EM will be performed whatever the value of argument \code{imputation}.
}
\note{Model indentification is not checked. Standard error and confidence intervals may be uncorrect if the model is not identified.}
\value{An object of class \code{dlsem}, with the following components:
  \item{estimate}{A list of objects of class \code{lm}, one for each response variable.}
  \item{model.code}{The model code after adaptation.}
  \item{exogenous}{The names of exogenous variables.}
  \item{group}{The name of the group factor. \code{NULL} is returned if \code{group}=\code{NULL}.}
  \item{log}{The value provided to argument \code{log}.}
  \item{ndiff}{The order of differentiation.}
  \item{data.orig}{The dataset provided to argument \code{data}.}
  \item{data.used}{Data used in the estimation, that is after eventual logarithmic transformation or differentiation.}
S3 methods available for class \code{dlsem} are:
  \item{print}{provides essential information on the structural model.}
  \item{summary}{shows summaries of estimation.}
  \item{plot}{displays the directed acyclic graph
  %(an edge is statistically significant if there is at least one time lag where the coefficient of
  %the variable originating the edge in the regression model of the variable receiving the edge is statistically significant)
  where each edge is coloured with respect to the sign of its causal effect (green: positive, red: negative, grey: not statistically significant).}
  \item{fitted}{returns fitted values.}
  \item{residuals}{returns residuals.}
  \item{predict}{returns predicted values.}
}
\references{A. Magrini, F. Bartolini, A. Coli, and B. Pacini (2016). Distributed-Lag Structural Equation Modelling:
An Application to Impact Assessment of Research Activity on European Agriculture.
\emph{Proceedings of the 48th Meeting of the Italian Statistical Society}, 8-10 June 2016, Salerno, IT.

W. K. Newey, and K. D. West (1978). A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix. \emph{Econometrica}, 55(3), 703-708.}
%\author{Alessandro Magrini <magrini@disia.unifi.it>}
\seealso{\link{unirootTest}, \link{applyDiff}, \link{EM.imputation}}
\examples{
data(industry)

# estimation without control options
mycode <- list(
  Consum~quec(Job,0,6),
  Pollution~quec(Job,1,11)+quec(Consum,1,6)
  )
myfit <- dlsem(mycode,group="Region",exogenous=c("Population","GDP"),
  data=industry,uniroot.check=TRUE,log=TRUE)


### adaptation of lag shapes (may take some seconds more)
## model code
#mycode <- list(
#  Consum~quec(Job,0,15),
#  Pollution~quec(Job,0,15)+quec(Consum,0,15)
#  )
#                      
#mycontrol <- list(
#  adapt=c(Consum=T,Pollution=T),
#  max.gestation=list(Consum=c(Job=3),Pollution=c(Consum=3,Job=3)),
#  min.width=list(Consum=c(Job=5),Pollution=c(Consum=5,Job=5)),
#  max.width=list(Consum=c(Job=15),Pollution=c(Consum=15,Job=15)),
#  sign=list(Consum=c(Job="+"),Pollution=c(Consum="+",Job="+"))
#  )
#myfit <- dlsem(mycode,group="Region",exogenous=c("Population","GDP"),data=industry,
#  control=mycontrol,uniroot.check=TRUE,log=TRUE)


# summaries of estimation
summary(myfit)

# display the directed acyclic graph
plot(myfit)
}
