\name{applyDiff}
\alias{applyDiff}
\title{Differentiation}
\description{Apply k-th order differentiation to one or more numerical variables. A single group factor may be taken into account.}
\usage{applyDiff(x = NULL, group = NULL, time = NULL, data, k = 0)}
\arguments{
  \item{x}{A vector containing the names of the numerical variables to be differenced. If \code{NULL}, all numerical variables will be differenced.}
  \item{group}{The name of a group factor (optional). If \code{NULL}, no groups are considered.}
  \item{time}{The name of a time factor (optional). If \code{NULL} and \code{group} is not \code{NULL}, data are assumed temporally ordered within each group.
  If both \code{group} and \code{time} are \code{NULL}, data are assumed temporally ordered.}
  \item{data}{An object of class \code{data.frame} containing the variables to be differenced,
  the group factor if \code{group} is not \code{NULL}, and the time factor if \code{time} is not \code{NULL}.}
  \item{k}{The order of differentiation to be applied to all variables, or a vector of orders of differentiation, one for each variable to be differenced.
  Default is 0 (no differentiation).}
  }
\value{The object of class \code{data.frame} provided to argument \code{data}, where the requested variables have been differenced.}
%\references{}
\author{Alessandro Magrini <magrini@disia.unifi.it>}
%\seealso{}
\examples{
data(agres)
applyDiff(c("NPATENT","GVA"),group="COUNTRY",data=agres,k=1)
}
