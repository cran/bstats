\name{white.test}
\alias{white.test}
\title{White test of constant variance}
\description{Perform a test to check the common variance assumption for a linear regression model.}
\usage{
  white.test(lmobj,squares.only=FALSE)
}
\arguments{
  \item{lmobj}{A fitted model}
  \item{squares.only}{If \code{TRUE}, include only the squared terms only.  
  Otherwise, include the interaction terms as well.}
}
\details{
  In this function, we check constant variance assmptions of the error terms.
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
 White, Halbert, A heteroskedasticity-consistent covariance 
 matrix estimator and a direct test for heteroskedasticity. 
 Econometrica 48 (1980), no. 4, 817-838.
}



\examples{

data(river)
lm0 = lm(Nitrogen~Agr+Forest+Rsdntial+ComIndl, data=river)
white.test(lm0)
white.test(lm0, squares.only=TRUE)

}
\keyword{htest}
