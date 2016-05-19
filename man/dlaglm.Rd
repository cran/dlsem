\name{dlaglm}
\alias{dlaglm}
\title{Distributed-lag linear regression}
\description{Fit a distributed-lag linear regression model. A single group factor may be taken into account using fixed effect panel estimation.}
\usage{dlaglm(formula, group = NULL, data, log = FALSE, L = 0, adapt = FALSE,
  max.gestation = NULL, min.width = NULL, sign = NULL)}
\arguments{
  \item{formula}{An object of class \code{formula} describing the model to be fitted. See \code{Details}.}
  \item{group}{The name of the group factor (optional). If \code{NULL}, no groups are considered.}
  \item{data}{An object of class \code{data.frame} containing the variables included in the model.}
  \item{log}{Logical. If \code{TRUE}, logarithmic transformation is applied to numerical variables. Default is \code{FALSE}.}
  \item{L}{The highest lag with non-zero autocorrelation (non-negative integer). If greater than 0, the Newey-West correction of the covariance matrix of estimates (Newey and West, 1987) is used. Default is 0.}
  \item{adapt}{Logical. If \code{TRUE}, AIC-based selection of lag shapes is performed (see \code{Details}). Default is \code{FALSE}.}
  \item{max.gestation}{See \code{Details}. Ignored if \code{adapt}=\code{FALSE}.}
  \item{min.width}{See \code{Details}. Ignored if \code{adapt}=\code{FALSE}.}
  \item{sign}{See \code{Details}. Ignored if \code{adapt}=\code{FALSE}.}
  }
\details{The formula cannot contain interaction terms (no ':' or '*' symbols), and
may contain the following operators for lag specification:

- \code{quec}: quadratic (2nd order polynomial) lag shape with endpoint constraints;

- \code{qdec}: quadratic (2nd order polynomial) decreasing lag shape.

Each operator must have the following three arguments (provided within brackets):

1) the name of the covariate to which the lag is applied;

2) the minimum lag with a non-zero coefficient;

3) the maximum lag with a non-zero coefficient.

For example, \code{quec(X1,3,15)} indicates that a quadratic lag shape with endpoint constraints must be applied to variable X1 in the interval (3,15).

If \code{adapt}=\code{TRUE}, AIC-based selection of lag shapes is performed.
The selection involves only the second and the third arguments of each lag operator in the model formula.
The following arguments of function \code{dlaglm} are considered in the selection:

- \code{max.gestation}: a named vector, including the maximum gestation lag for one or more covariates;

- \code{min.width}: a named vector, including the minimum lag width for one or more covariates;

- \code{sign}: a named vector, including the sign (either '+' or '-') of the coefficients of one or more covariates.
%
%For each of these arguments, values associated to covariates without a lag specification
%or to variables not included in the right side of the model formula will be ignored.
}
\value{An object of class \code{lm}.}
\note{If data contain unit roots, appropriate differentiation should be applied before.}
\references{
B. H. Baltagi (2008). Econometrics (Chapter 6). Springer Verlag, 4th edition, Berlin, DE.

W. K. Newey, and K. D. West (1978). A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix. \emph{Econometrica}, 55(3), 703-708.
}
\author{Alessandro Magrini <magrini@disia.unifi.it>}
\seealso{\link{unirootTest}, \link{applyDiff}, \link{EM.imputation}}
\examples{
data(agres)

# first order differentiation
agres.diff <- applyDiff(c("NPATENT","GVA","ENTR_INCOME","PPI"),group="COUNTRY",data=agres,k=1)

# fit without AIC-based selection of lag shapes
dl0 <- dlaglm(ENTR_INCOME~quec(NPATENT,0,3)+quec(GVA,1,10)+GDP+FARM_SIZE,group="COUNTRY",
  data=agres.diff,log=TRUE)

### AIC-based selection of lag shapes: NOT RUN
#dl0 <- dlaglm(ENTR_INCOME~quec(NPATENT,0,3)+quec(GVA,1,10)+GDP+FARM_SIZE,group="COUNTRY",
#  data=agres.diff,log=TRUE,adapt=T,max.gestation=c(NPATENT=3,GVA=3),
#  min.width=c(NPATENT=3,GVA=3),sign=c(NPATENT="+",GVA="+")) 

# summaries of model fitting
summary(dl0)
}
