\name{unirootTest}
\alias{unirootTest}
\title{Unit root test}
\description{Perform unit root test on a numerical variable. A single group factor may be taken into account.}
\usage{unirootTest(x, group = NULL, time = NULL, data, test = "adf", combine = "choi", k = 0)}
\arguments{
  \item{x}{The name of the numerical variable to be tested.}
  \item{group}{The name of a group factor (optional). If \code{NULL}, no groups are considered.}
  \item{time}{The name of a time factor (optional). If \code{NULL} and \code{group} is not \code{NULL}, data are assumed temporally ordered within each group.
  If both \code{group} and \code{time} are \code{NULL}, data are assumed temporally ordered.}
  \item{data}{An object of class \code{data.frame} containing the variables to be differenced,
  the group factor if \code{group} is not \code{NULL}, and the time factor if \code{time} is not \code{NULL}.}
  \item{test}{The test to use, that can be either \code{"adf"} (Augmented Dickey-Fuller, see Said and Dickey, 1984; Banerjee et al., 1993)
  or \code{"kpss"} (KPSS, see Kwiatkowski et al., 1992). Default is \code{"adf"}.}
  \item{combine}{The method to combine p-values of different groups, that can be either \code{"choi"} (Choi, 2001) or \code{"demetrescu"} (Demetrescu, 2006).
  Ignored if \code{group} is \code{NULL}. Default is \code{"choi"}.}
  \item{k}{The lag order to calculate the test statistic. Ignored if \code{test}=\code{"kpss"}. Default is 0.}
  }
\value{An object of class \code{htest}, containing the following components:
  \item{statistic}{The value of the test statistic.}
  \item{lag.order}{The lag order.}
  \item{alternative}{The alternative hypothesis.}
  \item{z.value}{The z-value of the test.}
  \item{p.value}{The p-value of the test.}
  \item{test}{The test used.}
  \item{combine}{The method used to combine p-values of different groups. \code{NULL} if \code{group} is \code{NULL}.}
  \item{n}{The number of observations per group.}
}
\note{The null hypothesis of the Augmented Dickey-Fuller test is the presence of an unit root. Instead, the null hypothesis of the KPSS test is the absence of an unit root.}
\references{
  A. Banerjee, J. J. Dolado, J. W. Galbraith, and D. F. Hendry (1993). Cointegration, Error Correction, and the Econometric Analysis of Non-Stationary Data. Oxford University Press, Oxford, UK.
  
  I. Choi (2001). Unit Root Tests for Panel Data. \emph{Journal of International Money and Finance}, 20, 249-272.
  
  M. Demetrescu, U. Hassler, and A. Tarcolea (2006). Combining Significance of Correlated Statistics with Application to Panel Data. \emph{Oxford Bulletin of Economics and Statistics}, 68(5), 647-663.
  
  D. Kwiatkowski, P. C. B. Phillips, P. Schmidt, Y. Shin (1992). Testing the Null Hypothesis of Stationarity against the Alternative of a Unit Root. \emph{Journal of Econometrics}, 54 (1-3): 159-178.
  
  S. E. Said, and D. A. Dickey (1984). Testing for Unit Roots in Autoregressive-Moving Average Models of Unknown Order. \emph{Biometrika}, 71, 599-607.
}
\author{Alessandro Magrini <magrini@disia.unifi.it>}
\seealso{\link{applyDiff}}
\examples{
data(agres)
unirootTest("NPATENT",group="COUNTRY",data=agres,k=1)
}
