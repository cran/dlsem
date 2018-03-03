\name{dlsem}
\alias{dlsem}
\title{Parameter estimation}
\description{Parameter estimation is performed for a distributed-lag linear structural equation model.
A single group factor may be taken into account.}
\usage{dlsem(model.code, group = NULL, time= NULL, exogenous = NULL, data, log = FALSE,
  diff.options = list(test="adf",combine="choi",k=0,lshort=TRUE,maxdiff=3),
  imput.options = list(tol=0.0001,maxiter=500,recons=TRUE,no.imput=NULL),
  global.control = NULL, local.control = NULL)}
\arguments{
  \item{model.code}{A list of objects of class \code{formula}, each describing a single regression model. See \code{Details}.}
  \item{group}{The name of the group factor (optional). If \code{NULL}, no groups are considered.}
  \item{time}{The name of the time factor (optional). If \code{time} is \code{NULL} and \code{group} is not \code{NULL}, data are assumed to be temporally ordered within each group.
  If both \code{time} and \code{group} are \code{NULL}, data are assumed to be temporally ordered.}
  \item{exogenous}{The name of exogenous variables (optional). Exogenous variables may be either quantitative or qualitative and must not appear in the model code.}  
  \item{data}{An object of class \code{data.frame} containing data.}
  \item{log}{Logical. If \code{TRUE}, logarithmic transformation is applied to strictly positive quantitative variables. Default is \code{FALSE}.}
  %
  \item{diff.options}{A list containing options for the differentiation. The list may consist of any number of components among the following: 
    \itemize{
    \item{test: }{the unit root test to use, that may be either \code{"adf"} or \code{"kpss"} (see \link{unirootTest}). Default is \code{"adf"};}
    \item{combine: }{the method to combine p-values of different groups, that may be either \code{"choi"} or \code{"demetrescu"} (see \link{unirootTest})
      Ignored if \code{group} is \code{NULL}. Default is \code{"choi"};}
    \item{k: }{the lag order to calculate the statistic of the Augmented Dickey-Fuller test.
      Ignored if \code{test}=\code{"kpss"}. Default is 0;}
    \item{lshort: }{logical. If \code{TRUE}, the short version of the truncation lag parameter is used for the KPSS test.
     Ignored if \code{test}=\code{"adf"}. Default is \code{TRUE};}
    \item{maxdiff: }{the maximum differentiation order to apply. If \code{maxdiff}=0, differentiation will not be applied. Default is 3.}
    }}
  %
  \item{imput.options}{A list containing options for the imputation of missing values through the Expectation-Maximization algorithm (Dempster et al., 1977).
  The list may consist of any number of components among the following:
    \itemize{
    \item{tol: }{the tolerance threshold. Default is 0.0001;}
    \item{maxiter: }{the maximum number of iterations. Default is 500. If \code{maxiter}=0, imputation will not be performed;}
    \item{recons: }{a logical value indicating whether univariate reconstruction of time series should be performed before applying differentiation. Reconstruction is performed using cubic splines and is limited to the observed range of each time series (no extrapolation). It may be useful to avoid further missing values due to differentiation in time series with contiguous missing values. Default is \code{TRUE}.} 
    \item{no.input: }{the name of variables to which reconstruction and imputation should not be applied.}  
    }}
  %
  \item{global.control}{A list containing global options for the estimation. The list may consist of any number of components among the following:
    \itemize{
    \item{adapt: }{a logical value indicating whether adaptation of lag shapes should be performed for all regression models. Default is \code{FALSE};}
    \item{max.gestation: }{the maximum gestation lag for all covariates. If not provided, it is assumed to be equal to \code{max.lead} (see below). Ignored if \code{adapt}=\code{FALSE};}
    \item{max.lead: }{the maximum lead lag for all covariates. If not provided, it is computed accordingly to the sample size. Ignored if \code{adapt}=\code{FALSE};} 
    \item{min.width: }{the minimum lag width for all covariates. It cannot be greater than \code{max.lead}. If not provided, it is assumed to be 0. Ignored if \code{adapt}=\code{FALSE};}
    \item{sign: }{the sign of parameter \eqn{\theta_j} (either '+' for positive or '-' for negative) for all covariates. If not provided, adaptation will disregard the sign of parameter \eqn{\theta_j}. Ignored if \code{adapt}=\code{FALSE};}  
    \item{selection: }{the criterion to be used for the adaptation of lag shapes, that may be one among \code{'aic'} (the default) to minimise the Akaike Information Criterion (Akaike, 1974), and
      \code{'bic'} to minimise the Bayesian Information Criterion (Schwarz, 1978).
      Default is \code{'aic'}.}
    }}
  
  \item{local.control}{A list containing variable-specific options for the estimation.
  These options prevail on the ones contained in \code{global.control}. See \code{Details}.}
}
\details{
Each regression model in a distributed-lag linear structural equation model has the form:
\deqn{y_t = \beta_0+\sum_{j=1}^J \sum_{l=0}^{L_j} \beta_{j,l} ~ x_{j,t-l}+\epsilon_t
%\hspace{1cm} %\epsilon_t \sim N(0,\sigma^2)
}
where \eqn{y_t} is the value of the response variable at time \eqn{t},
\eqn{x_{j,t-l}} is the value of the \eqn{j}-th covariate at \eqn{l} time lags before \eqn{t},
and \eqn{\epsilon_t} is the random error at time \eqn{t} uncorrelated with the covariates and with \eqn{\epsilon_{t-1}}.
The set \eqn{(\beta_{j,0},\beta_{j,1},\ldots,\beta_{j,L_j})}
is the lag shape of the \eqn{j}-th covariate.

For the endpoint-constrained quadratic lag shape:
\deqn{\beta_{j,l} = \theta_j \left[-\frac{4}{(b_j-a_j+2)^2} l^2+\frac{4(a_j+b_j)}{(b_j-a_j+2)^2} l-\frac{4(a_j-1)(b_j+1)}{(b_j-a_j+2)^2} \right] \hspace{1cm} a_j \leq l \leq b_j}
(otherwise, \eqn{\beta_{j,l}=0}).

For the quadratic decreasing lag shape:
\deqn{\beta_{j,l} = \theta_j \frac{l^2-2(b_j+1)l+(b_j+1)^2}{(b_j-a_j+1)^2} \hspace{1cm} a_j \leq l \leq b_j}
(otherwise, \eqn{\beta_{j,l}=0}).

For the gamma lag shape:
\deqn{\beta_{j,l} = \theta_j (l+1)^\frac{\delta}{1-\delta}\lambda_j^l \left[\left(\frac{\delta_j}{(\delta_j-1)\log(\lambda_j)}\right)^\frac{\delta_j}{1-\delta_j}\lambda_j^{\frac{\delta_j}{(\delta_j-1)\log(\lambda_j)}-1}\right]^{-1}}
\deqn{0<\delta_j<1 \hspace{1cm} 0<\lambda_j<1}

Formulas cannot contain qualitative variables or interaction terms (no ':' or '*' symbols), and
may contain the following operators for lag specification:
  \itemize{
  \item{\code{quec.lag}: }{quadratic (2nd order polynomial) lag shape with endpoint constraints;}
  \item{\code{qdec.lag}: }{quadratic (2nd order polynomial) decreasing lag shape;}
  \item{\code{gamma.lag}: }{gamma lag shape.}
  }
Each operator must have the following three arguments (provided within brackets):
  \enumerate{
  \item{the name of the covariate to which the lag is applied;}
  \item{parameter \eqn{a_j} (for 2nd order polynomial lag shapes), or parameter \eqn{\delta_j} (for the gamma lag shape);}
  \item{parameter \eqn{b_j} (for 2nd order polynomial lag shapes), or parameter \eqn{\lambda_j} (for the gamma lag shape).}
  }
For example, \code{quec.lag(X1,3,15)} indicates that a quadratic lag shape with endpoint constraints must be applied to variable X1 in the interval (3,15),
and \code{gamma.lag(X1,0.75,0.8)} indicates that a gamma lag shape with \eqn{\delta_j=0.75} and \eqn{\lambda_j=0.8} must be applied to variable X1.
See Judge et al. (1985, Chapters 9-10) for more details.

The formula of regression models with no endogenous covariates may be omitted from argument \code{model.code}.
The group factor and exogenous variables must not appear in any formula.

Argument \code{local.control} must be a named list containing one or more among the following components:
  \itemize{
  %\item{\code{L}: }{a named vector of non-negative integer values including the highest lag with non-zero autocorrelation for the indicated regression models.
  %If greater than 0, the Newey-West correction of the covariance matrix of estimates (Newey and West, 1987) is used.}
  \item{\code{adapt}: }{a named vector of logical values, where each component must have
the name of one endogenous variable and indicate whether adaptation of lag shapes should be
performed for the regression model of that variable.}
  \item{\code{max.gestation}: }{\code{max.gestation}: a named list. Each component of the list must have the
name of one endogenous variable and be a named vector.
Each component of the named vector must have the name of one covariate in the regression
model of the endogenous variable above and include the maximum gestation lag for its lag shape.}
  \item{\code{max.lead}: }{the same as \code{max.gestation}, with the exception that
the named vector must include the maximum lead lag.}
  \item{\code{min.width}: }{the same as \code{max.gestation}, with the exception that
the named vector must include the minimum lag width.}
  \item{\code{sign}: }{the same as \code{max.gestation}, with the exception that
the named vector must include the lag sign (either \code{'+'} for positive or \code{'-'} for negative).}
  Local control options have no default values, and global ones are applied in their absence.
  If some local control options conflict with global ones, only the former are applied.
  }
Imputation of missing values is performed after eventual logarithmic transformation
and differentiation, assuming group-specific means and time-invariant covariance matrix.

Qualitative variables cannot contain missing values.

Each quantitative variable must have at least a total of 3 observed values and,
if the group factor is specified, at least one observed value per group.
%
%Note that differentiation may create an additional number of missing values.
}
\value{An object of class \code{dlsem}, with the following components:
  \item{estimate}{A list of objects of class \code{lm}, one for each regression model.}
  \item{call}{A list containing the call for each regression model after eventual adaptation of lag shapes.}
  \item{exogenous}{The names of exogenous variables.}
  \item{group}{The name of the group factor. \code{NULL} is returned if \code{group}=\code{NULL}.}
  \item{log}{The value provided to argument \code{log}.}
  \item{uniroot}{the results of the uniroot test.}
  \item{ndiff}{The order of differentiation.}
  \item{diff.options}{Options used for the differentiation.}
  \item{imput.options}{Options used for the imputation of missing values.}
  \item{selection}{The criterion used for the adaptation of lag shapes.}
  \item{adaptation}{Variable-specific options used for the adaptation of lag shapes.}
  \item{Rsq}{The coefficient of determination (R-squared) of each regression model and the overall one.}  
  \item{data.orig}{The dataset provided to argument \code{data}.}
  \item{data.used}{Data used in the estimation, that is after eventual logarithmic transformation and differentiation.}
S3 methods available for class \code{dlsem} are:
  \item{print}{provides essential information on the model.}
  \item{summary}{shows summaries of estimation.}
  \item{plot}{displays the directed acyclic graph (DAG) of the model including only the endogenous variables.
  Option \code{conf} controls the confidence level (default is 0.95), while
  option \code{style} controls the style of the plot:
    \itemize{
    \item{\code{style=2} (the default): }{each edge is coloured with respect to the sign of the estimated causal effect (green: positive, red: negative, light grey: not statistically significant)};
    \item {\code{style=1}: }{edges with statistically significant causal effect are shown in black, otherwise they are shown in light grey};
    \item {\code{style=0}: }{all edges are shown in black disregarding statistical significance of causal effects}.
    }}
  \item{nobs}{returns the number of observations for each regression model.}
  \item{coef}{returns parameters \eqn{\theta_j (j=1,\ldots,J)} for each regression model.}
  \item{vcov}{returns the covariance matrix of estimates for each regression model.}
  \item{logLik}{returns the log-likelihood for each regression model.}
  \item{AIC}{returns the Akaike Information Criterion for each regression model and the overall one.}
  \item{BIC}{returns the Bayesian Information Criterion for each regression model and the overall one.}
  \item{fitted}{returns fitted values.}
  \item{residuals}{returns residuals.}
  \item{predict}{returns predicted values.
  If k-th order differentiation was applied, predicted values represent variations
  (on the original scale if \code{log}=\code{FALSE},
  percentage variations if \code{log}=\code{TRUE}) with respect to k periods before.
  Optionally, a data frame from which to predict may be provided to argument \code{newdata}.
  }
}
\references{
%W. K. Newey, and K. D. West (1978). A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix. \emph{Econometrica}, 55(3), 703-708.}

H. Akaike (1974). A New Look at the Statistical Identification Model. \emph{IEEE Transactions on Automatic Control}, 19, 716-723. DOI: 10.1109/TAC.1974.1100705

A. P. Dempster, N. M. Laird, and D. B. Rubin (1977). Maximum Likelihood from Incomplete Data via the EM Algorithm. \emph{Journal of the Royal Statistical Society}, Series B, 39(1): 1-38.

G. G. Judge, W. E. Griffiths, R. C. Hill, H. Lutkepohl, and T. C. Lee (1985). The Theory and Practice of Econometrics. John Wiley & Sons, 2nd ed., New York, US-NY.  ISBN: 978-0-471-89530-5

%P. Schmidt (1974). An Argument for the Usefulness of the Gamma Distributed Lag Model. \emph{International Economic Review}, 15(1).

G. Schwarz (1978). Estimating the Dimension of a Model. \emph{Annals of Statistics}, 6, 461-464. DOI: 10.1214/aos/1176344136
}
\seealso{\link{unirootTest}; \link{causalEff}.}
\examples{
data(industry)

# Estimation without adaptation of lag shapes
indus.code <- list(
  Consum~quec.lag(Job,0,4),
  Pollution~quec.lag(Job,2,6)+quec.lag(Consum,1,5)
  )
indus.mod <- dlsem(indus.code,group="Region",exogenous=c("Population","GDP"),data=industry,
  log=TRUE)

# Adaptation of lag shapes (estimation takes some seconds more)
indus.global <- list(adapt=TRUE,max.gestation=3,max.lead=15,min.width=5,sign="+")
indus.local <- list()
## NOT RUN:
# indus.mod <- dlsem(indus.code,group="Region",exogenous=c("Population","GDP"),data=industry,
#   global.control=indus.global,local.control=indus.local,log=TRUE)

# Summary of estimation
summary(indus.mod)

# DAG with coloured edges
plot(indus.mod)

# DAG disregarding statistical significance
plot(indus.mod,style=0)


### Comparison among alternative models

# Model 2: quadratic decreasing lag shapes
indus.code_2 <- list(
  Job ~ 1,
  Consum~qdec.lag(Job,0,15),
  Pollution~qdec.lag(Job,0,15)+qdec.lag(Consum,0,15)
  )
## NOT RUN:
# indus.mod_2 <- dlsem(indus.code_2,group="Region",exogenous=c("Population","GDP"),
#   data=industry,global.control=indus.global,local.control=indus.local,log=TRUE)

# Model 3: gamma lag shapes
indus.code_3 <- list(
  Job ~ 1,
  Consum~gamma.lag(Job,0.5,0.5),
  Pollution~gamma.lag(Job,0.5,0.5)+gamma.lag(Consum,0.5,0.5)
  )
## NOT RUN:
# indus.mod_3 <- dlsem(indus.code_3,group="Region",exogenous=c("Population","GDP"),
#   data=industry,global.control=indus.global,local.control=indus.local,log=TRUE)

# goodness of fit indices
## NOT RUN:
# lapply(list(QUEC=indus.mod,QDEC=indus.mod_2,GAMMA=indus.mod_3),AIC)
# lapply(list(QUEC=indus.mod,QDEC=indus.mod_2,GAMMA=indus.mod_3),BIC)


### A more complex model

data(agres)

# Qualitative exogenous variable
agres$POLICY <- factor(1*(agres$YEAR>=2005))
levels(agres$POLICY) <- c("no","yes")

# Causal levels
context.var <- c("GDP","EMPL_AGR","UAA","PATENT_OTHER","POLICY")
investment.var <- c("GBAORD_AGR","BERD_AGR")
research.var <- c("RD_EDU_AGR","PATENT_AGR")
impact.var <-  c("ENTR_INCOME_AGR","PPI_AGR")

# Constraints on lag shapes
agres.global <- list(adapt=TRUE,max.gestation=3,min.width=5,max.lead=20,sign="+")
auxcon1 <- rep(15,length(investment.var))
names(auxcon1) <- investment.var
auxcon2 <- rep("-",length(investment.var)+length(research.var))
names(auxcon2) <- c(investment.var,research.var)
agres.local <- list(max.lead=list(RD_EDU_AGR=auxcon1,PATENT_AGR=auxcon1),
  sign=list(PPI_AGR=auxcon2))

# Endpoint-constrained quadratic lag shapes (estimation takes a couple of minutes)
auxcode <- c(paste(investment.var,"~1",sep=""),
  paste(research.var,"~",paste("quec.lag(",investment.var,",0,20)",
    collapse="+",sep=""),sep=""),
  paste(impact.var,"~",paste("quec.lag(",c(investment.var,research.var),",0,20)",
    collapse="+",sep=""),sep=""))
agres.code <- list()
for(i in 1:length(auxcode)) {
  agres.code[[i]] <- formula(auxcode[i])
  }
## NOT RUN:
# agres.mod <- dlsem(agres.code,group="COUNTRY",exogenous=context.var,data=agres,
#   global.control=agres.global,local.control=agres.local,log=TRUE)
# summary(agres.mod)

## Gamma lag shapes (estimation takes a dozen of minutes)
auxcode_2 <- c(paste(investment.var,"~1",sep=""),
  paste(research.var,"~",paste("gamma.lag(",investment.var,",0.5,0.5)",
    collapse="+",sep=""),sep=""),
  paste(impact.var,"~",paste("gamma.lag(",c(investment.var,research.var),",0.5,0.5)",
    collapse="+",sep=""),sep=""))
agres.code_2 <- list()
for(i in 1:length(auxcode_2)) {
  agres.code_2[[i]] <- formula(auxcode_2[i])
  }
## NOT RUN:
# agres.mod_2 <- dlsem(agres.code_2,group="COUNTRY",exogenous=context.var,data=agres,
#  global.control=agres.global,local.control=agres.local,log=T)
# summary(agres.mod_2)
# lapply(list(QUEC=agres.mod,GAMMA=agres.mod_2),AIC)
# lapply(list(QUEC=agres.mod,GAMMA=agres.mod_2),BIC)
}
