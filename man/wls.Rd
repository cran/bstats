% This is a template R document file
% Distributed under GPL 3 or later

\name{wls}
\alias{wls}
\title{Weighted least squares estimate by groups}
\description{
Weighted least squares estimate by groups.
}
\usage{
wls(lmobj,group)
}
\arguments{
  \item{lmobj}{An R object by fitting an OLS model to a data set.}
  \item{group}{used to cluster the data.  Can be a factor or a numerical vector.}
}
\value{
  output the updated regressionn model with WLS.
}

\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}
\seealso{
 \code{\link{residual.plot}}.
}
\examples{

data(edu75)
lm0 = lm(Y~X1+X2+X3, data=edu75)
wls(lm0,group=edu75$Region)


}
\keyword{stats}

