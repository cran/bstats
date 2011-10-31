% This is a template R document file
% Distributed under GPL 3 or later

\name{lm.ci}
\alias{lm.ci}
\title{To compute the confidene interval of the regression parameters.}
\description{
To compute the confidene interval of the regression parameters.
}
\usage{
  lm.ci(lmobj,level=0.95)
}
\arguments{
  \item{lmobj}{An R object by fitting a linear regression model to a data set.}
  \item{level}{Confidence level. Default: 0.95.}
}

\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}

\seealso{
 \code{\link{model.test}}.
}

\examples{

data(ofc)
attach(ofc)
lm0 = lm(Head~Weight)
lm.ci(lm0)
lm1 = lm(Head~Weight+Gestation)
lm.ci(lm1, level=0.99)

}
\keyword{stats}

