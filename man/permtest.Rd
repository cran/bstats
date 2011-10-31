% This is a template R document file
% Distributed under GPL 3 or later

\name{perm.test}
\alias{perm.test}
\title{To perform a permutation test to compare two distributions.}
\description{
	To perform a permutation test to compare two distributions
}
\usage{
 perm.test(x, INDEX,iter=10000)
}
\arguments{
  \item{x}{Raw data.}
  \item{INDEX}{A vector indicating the group each value in X belong to.}
  \item{iter}{Iteration times, default: 10,000.}
}
\details{
To test a null hypothesis "H0: the two distributions are the same." against "H1: the two distributions are different".    
}
\value{

Missing values are not allowed.
}

\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}

\seealso{
 \code{\link{mediation.test}}.
}

\examples{

x1 = rnorm(100)
x2 = rnorm(200,3,9)
x = c(x1,x2)
group = c(rep(1,100),rep(2,200))
perm.test(x,group,iter=100)
}
\keyword{stats}

