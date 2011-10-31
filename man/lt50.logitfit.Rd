\name{ld50.logitfit}
\alias{ld50.logitfit}

\title{
Predict Doses for Binomial Assay model (using counts)
}
\description{
Calibrate binomial assays, generalizing the calculation of LD50 based on a logistic regression model.
}
\usage{
ld50.logitfit(rate, dose, p = 0.5)
}
\arguments{
\item{rate}{A vector of percentages of successes among all trials.}
\item{dose}{A vector of dosages.}
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
rate <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)/20

ld50.logitfit(rate,ldose,p = 0.5)
}
\keyword{regression}
\keyword{models}
