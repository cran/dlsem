\name{dlsem}
\alias{dlsem}
\title{Distributed-lag structural equation modelling}
\description{Estimation of a distributed-lag structural equation model with second-order polynomial lag shapes from data.}
\usage{dlsem(model.code, group = NULL, exogenous = NULL, data, log = FALSE, control = NULL,
  uniroot.check = TRUE, imputation = TRUE, test = "adf", combine = "choi", k = 0,
  lshort = TRUE, maxdiff = 5, tol = 0.0001, maxit = 500, plotDir = NULL)}
\arguments{
  \item{model.code}{A list of objects of class \code{formula}, each describing a single regression model. See \code{Details}.}
  \item{group}{The name of the group factor (optional). If \code{NULL}, no groups are considered.}
  \item{exogenous}{The name of exogenous variables (optional). Exogenous variables never appear on the left side of an equation and are not lagged.}  
  \item{data}{An object of class \code{data.frame} containing the variables included in the model.}
  \item{log}{Logical. If \code{TRUE}, logarithmic transformation is applied to quantitative variables. Default is \code{FALSE}.}
  \item{control}{A list containing options for estimation. See \code{Details}.}
  %\item{criterion}{The criterion for model selection, that can be either \code{"p.value"} to minimise the p-value, or \code{"aic"} to minimise the Akaike Information Criterion (AIC). Default is \code{"p.value"}.}
  \item{uniroot.check}{Logical. If \code{TRUE}, unit root test is performed for each variable, and appropriate differentation is applied. Default is \code{FALSE}.}
  \item{imputation}{Logical. If \code{TRUE}, missing data will be imputed using the EM algorithm. Default is \code{FALSE}.}
  \item{test}{The unit root test to use, that can be either \code{"adf"} or \code{"kpss"} (see \link{unirootTest}). Ignored if \code{uniroot.check}=\code{FALSE}. Default is \code{"adf"}.}
  \item{combine}{The method to combine p-values of different groups, that can be either \code{"choi"} or \code{"demetrescu"} (see \link{unirootTest}).
  Ignored if \code{uniroot.check}=\code{FALSE} or \code{group} is \code{NULL}. Default is \code{"choi"}.}
  \item{k}{The lag order to calculate the test statistic. Ignored if \code{test}=\code{"kpss"}. Default is 0.}
  \item{lshort}{Logical. If \code{TRUE}, the short version of the truncation lag parameter is used. Ignored if \code{test}=\code{"adf"}. Default is \code{TRUE}.}
  \item{maxdiff}{The maximum differentiation order to apply. Ignored if \code{uniroot.check}=\code{FALSE}. Default is 5.}
  \item{maxit}{The maximum number of iterations for the EM algorithm (see \link{EM.imputation}). Ignored if \code{imputation}=\code{FALSE}. Default is 500.}
  \item{tol}{The tolerance threshold of the EM algorithm (see \link{EM.imputation}). Ignored if \code{imputation}=\code{FALSE}. Default is 0.0001.}
  \item{plotDir}{A directory where to save the plots of the lag shapes (optional). If \code{NULL}, no plots will be produced.}
}
\details{Formulas cannot contain interaction terms (no ':' or '*' symbols), and
may contain the following operators for lag specification:

- \code{quec}: quadratic (2nd order polynomial) lag shape with endpoint constraints;

- \code{qdec}: quadratic (2nd order polynomial) decreasing lag shape.

Each operator must have the following three arguments (provided within brackets):

1) the name of the covariate to which the lag is applied;

2) the minimum lag with a non-zero coefficient;

3) the maximum lag with a non-zero coefficient.

For example, \code{quec(X1,3,15)} indicates that a quadratic lag shape with endpoint constraints must be applied to variable X1 in the interval (3,15).
The formula of regression models with no covariates excepting exogenous variables can be omitted from argument \code{model.code}.
Variables appearing in any formula are treated as quantitative.
The group factor and exogenous variables must appear in no formulas.

Argument \code{control} must be a named list containing one or more among the following components:

- \code{L}: a named vector of non-negative integer values including the highest lag with non-zero autocorrelation for one or more response variables.
  If greater than 0, the Newey-West correction of the covariance matrix of estimates (Newey and West, 1987) is used. Default is 0 for all response variables.

- \code{adapt}: a named vector of logical values indicating if adaptation of lag shapes must be performed for one or more response variables. Default is \code{FALSE} for all response variables.

- \code{max.gestation}: a named list. Each component of the list must refer to one response variable and contain a named vector, including the maximum gestation lag for one or more covariates.
  Ignored if \code{adapt}=\code{FALSE} for a certain covariate.

- \code{min.width}: a named list. Each component of the list must refer to one response variable and contain a named vector, including the minimum lag width for one or more covariates.
  Ignored if \code{adapt}=\code{FALSE} for a certain covariate.

- \code{sign}: a named list. Each component of the list must refer to one response variable and contain a named vector, including the sign
  (either '+' for non-negative, or '-' for non-positive) of the coefficients of one or more covariates. Ignored if \code{adapt}=\code{FALSE} for a certain covariate.
%
%For the last three components, values associated to covariates without a lag specification
%or to variables not appearing in the model code will be ignored.

Variables appearing in the model code but not included in the dataset will be considered as unobserved.
If there is at least one unobserved variable, imputation using EM will be performed whatever the value of argument \code{imputation}.
}
\note{Model indentification is not checked. Standard error and confidence intervals may be uncorrect if the model is not identified.}
\value{An object of class \code{dlsem}, with the following components:
  \item{estimate}{A list of objects of class \code{lm}, one for each response variable.}
  \item{model.code}{The model code.}
  \item{exogenous}{The names of exogenous variables.}
  \item{group}{The name of the group factor. \code{NULL} is returned if \code{group}=\code{NULL}.}
  \item{data}{The data used for estimation.}
S3 methods available for class \code{dlsem} are:
  \item{print}{provides essential information on the structural model.}
  \item{summary}{shows summaries of estimation.}
  \item{plot}{displays the directed acyclic graph including only statistically significant edges
  (an edge is statistically significant if there is at least one time lag where the coefficient of
  the variable originating the edge in the regression model of the variable receiving the edge
  is statistically significant). The method contains a further argument:
  \code{sign.col}, a logical value indicating whether edges should be coloured with respect to their sign (green: positive, red: negative. Default is \code{TRUE}).}
  \item{fitted}{returns fitted values.}
  \item{residuals}{returns residuals.}
  \item{predict}{returns predicted values.}
  \item{deviance}{returns the deviance of the structural model.}
  \item{extractAIC}{computes the Akaike Information Criterion (AIC) of the structural model.}
}
\references{A. Magrini, F. Bartolini, A. Coli, and B. Pacini (2016). Distributed-Lag Structural Equation Modelling:
An Application to Impact Assessment of Research Activity on European Agriculture.
\emph{Proceedings of the 48th Meeting of the Italian Statistical Society}, 8-10 June 2016, Salerno, IT.

W. K. Newey, and K. D. West (1978). A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix. \emph{Econometrica}, 55(3), 703-708.}
\author{Alessandro Magrini <magrini@disia.unifi.it>}
\seealso{\link{dlsem-package}, \link{unirootTest}, \link{applyDiff}, \link{EM.imputation}}
\examples{
data(agres)

# estimation without control options
mycode <- list(
  GVA~quec(NPATENT,1,15),
  PPI~quec(NPATENT,0,13)+quec(GVA,0,14),
  ENTR_INCOME~quec(NPATENT,0,14)+quec(GVA,1,14)
  )
myfit <- dlsem(mycode,group="COUNTRY",exogenous=c("GDP","FARM_SIZE"),data=agres,
  uniroot.check=TRUE,imputation=FALSE,log=TRUE)


### adaptation of lag shapes (may take some seconds more)
## model code
#mycode <- list(
#  GVA~quec(NPATENT,0,15),
#  PPI~quec(NPATENT,0,15)+quec(GVA,0,15),
#  ENTR_INCOME~quec(NPATENT,0,15)+quec(GVA,0,15)
#  )
#
## control options
#mycontrol <- list(
#  adapt=c(GVA=TRUE,PPI=TRUE,ENTR_INCOME=TRUE),
#  max.gestation=list(GVA=c(NPATENT=3),PPI=c(NPATENT=3,GVA=3),ENTR_INCOME=c(NPATENT=3,GVA=3)),
#  min.width=list(GVA=c(NPATENT=5),PPI=c(NPATENT=5,GVA=5),ENTR_INCOME=c(NPATENT=5,GVA=5)),
#  sign=list(GVA=c(NPATENT="+"),PPI=c(NPATENT="-",GVA="-"),ENTR_INCOME=c(NPATENT="+",GVA="+"))
#  )
#
#myfit <- dlsem(mycode,group="COUNTRY",exogenous=c("GDP","FARM_SIZE"),data=agres,
#  control=mycontrol,uniroot.check=TRUE,imputation=TRUE,log=TRUE)


# summaries of estimation
summary(myfit)

# display the DAG with significant edges only
plot(myfit)
}
