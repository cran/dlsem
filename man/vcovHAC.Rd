\name{vcovHAC}
\alias{vcovHAC}
\title{Heteroskedasticty and autocorrelation consistent covariance matrix}
\description{The heteroskedasticty and autocorrelation consistent (HAC) covariance matrix of least square estimates (Newey & West, 1978) is computed. A single group factor may be taken into account.}
\usage{vcovHAC(x, group = NULL)}
\arguments{
  \item{x}{An object of class \code{lm} or \code{dlsem}.}
  \item{group}{The name of the group factor (optional). If \code{NULL}, no groups are considered.}
  }
\value{A matrix if \code{x} is of class \code{lm}, or, if \code{x} is of class \code{dlsem}, a list of matrices, one for each regression model.
Each matrix has the attribute \code{max.lag}, indicating the maximum lag of autocorrelation, automatically computed based on fit to data.
If \code{group} is not \code{NULL}, this is computed within each group.}
\note{Residuals are assumed to be temporally ordered within each group.}
\references{
W. K. Newey, and K. D. West (1978). A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix. \emph{Econometrica}, 55(3), 703-708.
}
%\seealso{}
\examples{
data(industry)
m0 <- lm(Consum ~ -1+Region+quec.lag(Job,0,5,x.group=Region), data=industry)
# traditional covariance matrix
S <- vcov(m0)
# HAC covariance matrix
S_hac <- vcovHAC(m0,group="Region")
# comparison of standard errors
cbind(traditional=sqrt(diag(S)),hac=sqrt(diag(S_hac)))
}
