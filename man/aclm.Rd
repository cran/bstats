\name{ac}
\alias{ac}
\alias{ac.lm}
\alias{ac.default}

\title{Autocorrelation}
\description{
  Removal of autocorrelation by transformation.
}
\usage{
ac(lmobj,type='cochrane', ...)

\method{ac}{lm}(lmobj,type='cochrane', ...)
}
 
\arguments{
  \item{lmobj}{an object that inherits from class \code{lm}, such as an \code{lm}
  or \code{glm} object.}
  \item{type}{method selection: 'iterative', 'cochrane'.}
  \item{\dots}{not used.}
}
\details{
	'iterative': simultaneously estimate the regression coefficients and rho by minimizing the sum squared errors. A grid searching method is used.

'cochrane': 1. Fit a linear regression model and compute OLS estimates 2. Calculate the residuals to estimate rho from the data. 3. Fit (1) to obtain estimates of the regression coefficients.
4. Check to see whether autocorrelation still exist. If yes, repeat by using the estimated coefficients from step 3 in step 1.
}
\value{
	coefficients, rhohat, dwtest, re-fitted model.
}
\references{ 
Cochrane and Orcutt (1949)

St 335 text
}
\author{Wang, B.}

\examples{

data(edu75)
lm0 = lm(Y~X1+X2+X3, data=edu75)
ac.lm(lm0,type='iterative')
ac.lm(lm0, type='cochrane')

}

\keyword{regression}
