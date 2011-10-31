% This is a template R document file
% Distributed under GPL 3 or later

\name{edf}
\alias{edf}
\alias{print.edf}
\alias{plot.edf}
\title{To compute the empirical distribution function.}
\description{
To compute the empirical distribution function.}
\usage{
  edf(x,y=NULL)
}
\arguments{
  \item{x}{A sample. 'NA' values will be automatically removed.}
  \item{y}{A grid of points where the edf will be evaluated.}
}
\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}
\seealso{
 \code{\link{scb}}.
}
\examples{

x = rnorm(100)
(out = edf(x))
plot(out)
(out2= scb(out))
lines(out2)

}
\keyword{stats}

