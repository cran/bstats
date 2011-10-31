% This is a template R document file
% Distributed under GPL 3 or later

\name{residual.plot}
\alias{residual.plot}
\title{Draw residual plots for an ordinary regression model.}
\description{
Draw residual plots for an ordinary regression model.
}
\usage{
  residual.plot(lmobj,type='fitted',col=1)
}
\arguments{
  \item{lmobj}{An R object by fitting an OLS model to a data set.}
  \item{type}{Type of residual plot(s): 
  'fitted', residuals against fitted values; 
  'index', residuals against index; 
  'predictor', residuals against each of the predictors in the fitted model; 
  'qqplot', qq-plot of the standardized residuals to check the normality assumption.}
  \item{col}{Color of the plot.}
}
\value{
Missing values are not allowed.
}

\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}
\seealso{
 \code{\link{influential.plot}}.
}
\examples{

data(river)
lm0 = lm(Nitrogen~Agr+Forest+Rsdntial+ComIndl, data=river)
residual.plot(lm0)
residual.plot(lm0,type='index')
residual.plot(lm0,type='predictor')


}
\keyword{stats}

