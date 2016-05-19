\name{dlsem}
\alias{dlsem}
\title{Distributed-Lag Structural Equation Modelling.}
\description{Fit a distributed-lag structural equation model.}
\usage{dlsem(model.code, group = NULL, context = NULL, data, log = FALSE, control = NULL,
  uniroot.check = TRUE, imputation = TRUE, test = "adf", combine= "choi", maxdiff = 5,
  tol = 0.0001, maxit = 500, plotDir = NULL)}
\arguments{
  \item{model.code}{A list of objects of class \code{formula}, each describing a single regression model. See \code{Details}.}
  \item{group}{The name of the group factor (optional). If \code{NULL}, no groups are considered.}
  \item{context}{The name of context variables (optional). Context variables never appear on the right side of an equation and are not lagged.}  
  \item{data}{An object of class \code{data.frame} containing the variables included in the model.}
  \item{log}{Logical. If \code{TRUE}, logarithmic transformation is applied to numerical variables. Default is \code{FALSE}.}
  \item{control}{A list containing various options for the fitting. See \code{Details}.}
  \item{uniroot.check}{Logical. If \code{TRUE}, unit root test is performed for each variable, and appropriate differentation is applied. Default is \code{FALSE}.}
  \item{imputation}{Logical. If \code{TRUE}, missing data will be imputed using the EM algorithm. Default is \code{FALSE}.}
  \item{test}{The unit root test to use, that can be either \code{"adf"} or \code{"kpss"} (see \link{unirootTest}). Ignored if \code{uniroot.check}=\code{FALSE}. Default is \code{"adf"}.}
  \item{combine}{The method to combine p-values of different groups, that can be either \code{"choi"} or \code{"demetrescu"} (see \link{unirootTest}).
  Ignored if \code{uniroot.check}=\code{FALSE} or \code{group} is \code{NULL}. Default is \code{"choi"}.}
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
The formula of regression models with no covariates excepting context variables can be omitted from argument \code{model.code}.

Argument \code{control} must be a named list containing one or more among the following components:

- \code{L}: a named vector of non-negative integer values including the highest lag with non-zero autocorrelation for a certain response variable. Default is 0 for all response variables.

- \code{adapt}: a named vector of logical values indicating if AIC-based selection of lag shapes must be performed for a certain response variable. Default is \code{FALSE} for all response variables.

- \code{max.gestation}: a named list. Each component of the list must refer to one response variable and contain a named vector, including the maximum gestation lag for one or more covariates.

- \code{min.width}: a named list. Each component of the list must refer to one response variable and contain a named vector, including the minimum lag width for one or more covariates.

- \code{sign}: a named list. Each component of the list must refer to one response variable and contain a named vector, including the sign (either '+' or '-') of the coefficients of one or more covariates.
%
%For the last three components, values associated to covariates without a lag specification
%or to variables not mentioned in the model code will be ignored.
}
\value{An object of class \code{dlsem}, with the following components:
  \item{estimate}{A list of objects of class \code{lm}, one for each response variable.}
  \item{model.code}{The model code.}
  \item{context}{The names of context variables.}
  \item{group}{The name of the group factor. \code{NULL} is returned if \code{group}=\code{NULL}.}
S3 methods available for class \code{dlsem} are:
  \item{print}{provides essential information on the structural model.}
  \item{summary}{shows summaries of model fitting.}
  \item{plot}{displays the directed acyclic graph of the structural model. The method contains a further argument:
  \code{sign.col}, a logical value indicating whether edges should be coloured with respect to their sign (green: positive, red: negative. Default is \code{TRUE}).}
  \item{fitted}{returns fitted values.}
  \item{residuals}{returns residuals.}
  \item{predict}{returns predicted values.}
  \item{extractAIC}{computes the Akaike Information Criterion (AIC) of the structural model.}
}
\references{A. Magrini, F. Bartolini, A. Coli, and B. Pacini, 2016. Distributed-Lag Structural Equation Modelling:
An Application to Impact Assessment of Research Activity on European Agriculture.
To be appeared in: \emph{Proceedings of the 48th Meeting of the Italian Statistical Society}, 8-10 June 2016, Salerno, IT.}
\author{Alessandro Magrini <magrini@disia.unifi.it>}
\seealso{\link{dlsem-package}, \link{unirootTest}, \link{applyDiff}, \link{EM.imputation}}
\examples{
data(agres)

# model code
mycode <- list(
  GVA~quec(NPATENT,0,4),
  PPI~quec(NPATENT,0,6)+quec(GVA,0,10),
  ENTR_INCOME~quec(NPATENT,0,3)+quec(GVA,1,10)
  )

# control options
mycontrol <- list(
  adapt=c(GVA=TRUE,PPI=TRUE,ENTR_INCOME=TRUE),
  max.gestation=list(GVA=c(NPATENT=3),PPI=c(NPATENT=3,GVA=3),ENTR_INCOME=c(NPATENT=3,GVA=3)),
  min.width=list(GVA=c(NPATENT=3),PPI=c(NPATENT=3,GVA=3),ENTR_INCOME=c(NPATENT=3,GVA=3)),
  sign=list(GVA=c(NPATENT="+"),PPI=c(NPATENT="-",GVA="-"),ENTR_INCOME=c(NPATENT="+",GVA="+"))
  )

# fitting without control options
mod0 <- dlsem(mycode,group="COUNTRY",context=c("GDP","FARM_SIZE"),data=agres,
  uniroot.check=TRUE,maxdiff=1,imputation=FALSE,log=TRUE)

### AIC-based selection of lag shapes: not RUN
#mod0 <- dlsem(mycode,group="COUNTRY",context=c("GDP","FARM_SIZE"),data=agres,
#  control=mycontrol,uniroot.check=TRUE,maxdiff=1,imputation=FALSE,log=TRUE)

# summaries of model fitting
summary(mod0)

# the directed acyclic graph of the structural model
plot(mod0)
}
