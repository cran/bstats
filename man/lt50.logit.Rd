\name{ld50.logit}
\alias{ld50.logit}
\alias{print.glm.dose}
\alias{plot.glm.dose}
\alias{lines.glm.dose}

\title{
Predict Doses for Binomial Assay model (using counts)
}
\description{
Calibrate binomial assays, generalizing the calculation of LD50 based on a logistic regression model.
}
\usage{
ld50.logit(ndead, ntotal, dose, cf = 1:2, p = 0.5)
}
\arguments{
\item{ndead}{A vector of number of failures.}
\item{ntotal}{Total number of trials.}
\item{dose}{A vector of dosages.}
\item{cf}{
The terms in the coefficient vector giving the intercept and
coefficient of (log-)dose
}
\item{p}{
Probabilities at which to predict the dose needed.
}}
\references{
  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S.}
  Springer.
}
\examples{
ldose <- rep(0:5, 2)
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
n=20

ld50.logit(numdead,n,ldose,p = 0.5)
}
\keyword{regression}
\keyword{models}
