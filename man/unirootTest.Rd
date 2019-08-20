\name{unirootTest}
\alias{unirootTest}
\title{Unit root test}
\description{The Augmented Dickey-Fuller test is performed on a set of quantitative variables. A single group factor may be taken into account.}
\usage{unirootTest(x = NULL, group = NULL, time = NULL, data, combine = "choi", log = FALSE)}
\arguments{
  \item{x}{A vector including the name of the quantitative variables to be tested.
  If \code{NULL} (the default), all the quantitative variables contained in \code{data} will be tested.}
  \item{group}{The name of the group factor (optional). If \code{NULL}, no groups are considered.}
  \item{time}{The name of the time factor (optional).  This variable must be either a numeric identificative or a date in format '\%Y/\%m/\%d','\%d/\%m/\%Y', or '\%Y-\%m-\%d'. If \code{time} is \code{NULL} and \code{group} is not \code{NULL}, data are assumed to be temporally ordered within each group.
  If both \code{time} and \code{group} are \code{NULL}, data are assumed to be temporally ordered.}
  \item{data}{An object of class \code{data.frame} containing the variables to be tested,
  the group factor if \code{group} is not \code{NULL}, and the time factor if \code{time} is not \code{NULL}.}
  \item{combine}{The method to combine p-values of different groups, that may be either \code{"choi"} (Choi, 2001) or \code{"demetrescu"} (Demetrescu, 2006).
  Ignored if \code{group} is \code{NULL}. Default is \code{"choi"}.}
  \item{log}{Logical. If \code{TRUE}, logarithmic transformation is applied to all strictly positive quantitative variables. Default is \code{FALSE}.}
  }
\value{An object of class \code{unirootTest}, consisting of a list with one component for each variable tested. Each list contains the following components:
  \item{statistic}{The value of the test statistic.}
  \item{lag.order}{The lag order at which the test statistic is computed. It is automatically selected according to the precedure by Ng \& Perron (2001).}
  \item{n}{The total number of observations if \code{group} is \code{NULL}, otherwise the number of observations per group.}
  \item{z.value}{The z-value of the test.}
  \item{p.value}{The p-value of the test.}
}
\note{The null hypothesis of the Augmented Dickey-Fuller test is the presence of a unit root.
The lag order to calculate the test statistic is automatically selected according to the precedure by Ng \& Perron (2001).}
\references{
I. Choi (2001). Unit Root Tests for Panel Data. \emph{Journal of International Money and Finance}, 20, 249-272. DOI: 10.1016/S0261-5606(00)00048-6

M. Demetrescu, U. Hassler, and A. Tarcolea (2006). Combining Significance of Correlated Statistics with Application to Panel Data. \emph{Oxford Bulletin of Economics and Statistics}, 68(5), 647-663. DOI: 10.1111/j.1468-0084.2006.00181.x
  
D. A. Dickey, and W. A. Fuller (1981). Likelihood Ratio Statistics for Autoregressive Time Series with a Unit Root. \emph{Econometrica}, 49: 1057-1072. DOI: 10.2307/1912517

S. Ng, and W. P. Perron (2001). Lag Length  Selection and the Construction of Unit Root Tests with Good Size and Power. \emph{Econometrica}, 60: 1519-1554. DOI: 10.1111/1468-0262.00256.
}
%\seealso{}
\examples{
data(industry)
indus.urt <- unirootTest(c("Job","Consum","Population","GDP"),
  group="Region",time="Year",data=industry)
indus.urt      ## p-values
indus.urt$Job  ## details for variable 'Job'
}
