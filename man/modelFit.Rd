\name{modelFit}
\alias{modelFit}
\title{Fit indices}
\description{Compute local and global fit indices.}
\usage{modelFit(x)}
\arguments{
  \item{x}{An object of class \code{dlsem}.}
  }
\value{A list of five numeric vectors.
The fist component contains the R squared for each regression model and the overall one.
The second component contains the log-likelihood for each regression model and the overall one.
The third component contains the Akaike Information Criterion (Akaike, 1974) for each regression model and the overall one.
The fourth component contains the Bayesian Information Criterion (Schwarz, 1978) for each regression model and the overall one.
The fifth component contains the Minimum Information Length (Rissanen, 1978) for each regression model and the overall one.
}
%\note{}
\references{
H. Akaike (1974). A New Look at the Statistical Identification Model. \emph{IEEE Transactions on Automatic Control}, 19, 716-723. DOI: 10.1109/TAC.1974.1100705

J. Rissanen (1978). Modeling by Shortest Data Description. \emph{Automatica}, 14(5): 465-658. DOI: 10.1016/0005-1098(78)90005-5

G. Schwarz (1978). Estimating the Dimension of a Model. \emph{Annals of Statistics}, 6, 461-464. DOI: 10.1214/aos/1176344136
}
\seealso{\link{dlsem}}
\examples{
data(industry)
mycode <- list(
  Consum~quec(Job,0,5),
  Pollution~quec(Job,1,8)+quec(Consum,1,6)
  )
myfit <- dlsem(mycode,group="Region",exogenous=c("Population","GDP"),data=industry,log=TRUE)
modelFit(myfit)
}
