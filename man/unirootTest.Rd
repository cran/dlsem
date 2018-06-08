\name{unirootTest}
\alias{unirootTest}
\title{Unit root test}
\description{The Augmented Dickey-Fuller test is performed on a set of quantitative variables. A single group factor may be taken into account.}
\usage{unirootTest(x, group = NULL, time = NULL, data, combine = "choi", k = 0)}
\arguments{
  \item{x}{A vector including the name of the quantitative variables to be tested.}
  \item{group}{The name of the group factor (optional). If \code{NULL}, no groups are considered.}
  \item{time}{The name of the time factor (optional). If \code{time} is \code{NULL} and \code{group} is not \code{NULL}, data are assumed to be temporally ordered within each group.
  If both \code{time} and \code{group} are \code{NULL}, data are assumed to be temporally ordered.}
  \item{data}{An object of class \code{data.frame} containing the variables to be differenced,
  the group factor if \code{group} is not \code{NULL}, and the time factor if \code{time} is not \code{NULL}.}
  \item{combine}{The method to combine p-values of different groups, that may be either \code{"choi"} (Choi, 2001) or \code{"demetrescu"} (Demetrescu, 2006).
  Ignored if \code{group} is \code{NULL}. Default is \code{"choi"}.}
  \item{k}{The lag order to calculate the statistic of the Augmented Dickey-Fuller test. If \code{k=0} (the default), the standard Dickey-Fuller test is performed.
  If \code{k=NULL}, the lag order is taken as the cubic root of the sample size.}
  }
\value{An object of class \code{unirootTest}, consisting of a list with one component for each quantitative variable. Each list contains the following components:
  \item{statistic}{The value of the test statistic.}
  \item{n}{The total number of observations if \code{group} is \code{NULL}, otherwise the number of observations per group.}
  \item{z.value}{The z-value of the test.}
  \item{p.value}{The p-value of the test.}
}
\note{The null hypothesis of the Augmented Dickey-Fuller test is the presence of a unit root.}
\references{
I. Choi (2001). Unit Root Tests for Panel Data. \emph{Journal of International Money and Finance}, 20, 249-272. DOI: 10.1016/S0261-5606(00)00048-6

M. Demetrescu, U. Hassler, and A. Tarcolea (2006). Combining Significance of Correlated Statistics with Application to Panel Data. \emph{Oxford Bulletin of Economics and Statistics}, 68(5), 647-663. DOI: 10.1111/j.1468-0084.2006.00181.x
  
D. A. Dickey, and W. A. Fuller (1981). Likelihood Ratio Statistics for Autoregressive Time Series with a Unit Root. \emph{Econometrica}, 49: 1057-1072. DOI: 10.2307/1912517
%    
%D. Kwiatkowski, P. C. B. Phillips, P. Schmidt, and Y. Shin (1992). Testing the Null Hypothesis of Stationarity against the Alternative of a Unit Root. \emph{Journal of Econometrics}, 54 (1-3): 159-178. DOI: 10.1016/0304-4076(92)90104-Y
}
%\seealso{}
\examples{
data(industry)
unirootTest(c("Job","Consum","GDP"),group="Region",data=industry)
}
