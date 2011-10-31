% This is a template R document file
% Distributed under GPL 3 or later

\name{scb}
\alias{scb}
\alias{plot.scb}
\alias{print.scb}
\alias{lines.scb}
\title{To compute the simultaneous confidence bands.}
\description{
To compute the simultaneous confidence bands.}
\usage{
  scb(x,alpha=0.05)
}
\arguments{
  \item{x}{An R object.  Currently, only 'edf' objects are supported.}
  \item{alpha}{Significance level.  Default 0.05 for a 95 percent confidence level.}
}
\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}
\seealso{
 \code{\link{edf}}.
}
\examples{

x = rnorm(100)
(out = edf(x))
plot(out)
(out2= scb(out))
lines(out2)

}
\keyword{stats}

