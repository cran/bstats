% This is a template R document file
% Distributed under GPL 3 or later

\name{model.test}
\alias{model.test}
\title{To compare two models and determine which one is adequate.}
\description{
To compare a full model and reduced model to test whether the reduced model is adequate or not.
}
\usage{
  model.test(fmobj,rmobj,alpha=0.05)
}
\arguments{
  \item{fmobj}{An R object by fitting a full linear regression model (FM) to a data set.}
  \item{rmobj}{An R object by fitting a reduced linear regression model (RM) to a data set.}
  \item{alpha}{Significance level. Default: alpha=0.05.}
}
\details{
To test a null hypothesis "H0: the RM is adequate" against "H1: the FM is adequate".  The values of test statistic, p-value and critical value based on an F test will be given.  
}
\value{

Missing values are not allowed.
}

\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\seealso{
 \code{\link{lm.ci}}.
}

\examples{

data(supervisor)
attach(supervisor)
lm0 = lm(Y~X1+X3)
lm1 = lm(Y~X1+X2+X3+X4+X5+X6)
model.test(lm1,lm0)

}
\keyword{stats}

