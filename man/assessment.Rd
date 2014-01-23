% This is a template R document file
% Distributed under GPL 3 or later

\name{assessment}
\alias{assessment}

\title{Assessment Summaries}
\description{
 To compute summary statistics for assessment.}
\usage{
  assessment(x, keys, cutoff=0.60)
}
\arguments{
  \item{x}{A data frame of the answers for the common problems for assessment.}
  \item{keys}{Solution keys for the common problems.}
  \item{cutoff}{A cutoff to determine whether a student pass the criteria (\code{>= cutoff}) 
  or not (\code{<cutoff}).}
}
\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\examples{

  data(st210107)
  keys <- c("C", "A", "B", "C", "D")
  out <- assessment(st210107, keys)
  out
}
\keyword{stats}

