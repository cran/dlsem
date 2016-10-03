\name{EM.imputation}
\alias{EM.imputation}
\title{Imputation of missing values using Expectation-Maximization}
\description{Missing values of quantitative variables are replaced by their expected value computed using the Expectation-Maximization (EM) algorithm (Dempster et al., 1977). A Multivariate Gaussian distribution is assumed.}
\usage{EM.imputation(x = NULL, group = NULL, data, plotDir = NULL, tol = 0.0001, maxit = 500)}
\arguments{
  \item{x}{A vector containing the names of the quantitative variables to be imputed. If \code{NULL}, all quantitative variables will be imputed.}
  \item{group}{The name of a group factor (optional). If \code{NULL}, no groups are considered.}
  \item{data}{An object of class \code{data.frame} containing the variables to be imputed and, if \code{group} is not \code{NULL}, the group factor.}
  \item{plotDir}{A directory where to save the plots after imputation (optional). If \code{NULL}, no plots will be produced.}
  \item{tol}{The tolerance threshold of the EM algorithm. Default is 0.0001.}
  \item{maxit}{The maximum number of iterations for the EM algorithm. Default is 500.}
}
\value{The object of class \code{data.frame} provided to argument \code{data}, where the requested variables have been imputed.}
\note{If data contain unit roots, appropriate differentiation should be applied before.}
\references{A. P. Dempster, N. M. Laird, D. B. Rubin (1977). Maximum Likelihood from Incomplete Data via the EM Algorithm. \emph{Jornal of the Royal Statustical Society}, Series B, 39(1), 1-38.}
%\author{Alessandro Magrini <magrini@disia.unifi.it>}
\seealso{\link{unirootTest}, \link{applyDiff}}
\examples{
data(agres)
EM.imputation(c("PPI","ENTR_INCOME"),group="COUNTRY",data=agres)
}
