% This is a template R document file
% Distributed under GPL 3 or later

\name{eda}
\alias{eda}
\title{Exploratory Data Analysis (EDA)}
\description{
A graphical method for exploratory data analysis.
}
\usage{
  eda(x)
}
\arguments{
  \item{x}{A vector of numeric values. 'NA' values will be automatically removed.}
}
\author{
  B. Wang  \email{bwang@southalabama.edu}
}
\details{
  This function is modified based on an S-Plus function written by Dr. Jiayang Sun 
  (Course ST 426. Statistics Department, Case Western Reserve University).
.
}
\references{
 J. W. Tukey. (1977). Exploratory Data Analysis. Addison-Wesley, New York.
}
\seealso{
 \code{\link{scb}}.
}
\examples{
 x = rnorm(100)
 eda(x)
}
\keyword{stats}

