\name{model.check}
\alias{model.check}
\title{Linear Regression Model Check}
\description{Performs tests to check the least squares assumptions for a linear regression model.}
\usage{
  model.check(lmobj)
}
\arguments{
  \item{lmobj}{A fitted model}
}
\details{
  In this function, we check the normality, independece, and constant variance assmptions of the error terms, and the presence of multicollinearity.
}

\value{
 A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{p.value}{the p-value of the test.}
  \item{parameter}{degrees of freedom.}
  \item{method}{a character string indicating what type of test was
    performed.}
  \item{data.name}{a character string giving the name(s) of the data.}
}

\references{
To be updated.
}


\examples{

data(river)
lm0 = lm(Nitrogen~Agr+Forest+Rsdntial+ComIndl, data=river)
model.check(lm0)

}
\keyword{htest}
