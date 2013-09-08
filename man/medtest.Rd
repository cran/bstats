% This is a template R document file
% Distributed under GPL 3 or later

\name{mediation.test}
\alias{mediation.test}
\title{The Sobel mediation test}
\description{
  To compute statistics and p-values for the Sobel test.  Results for three versions of "Sobel test" are provided: Sobel test, Aroian test and Goodman test. 
}
\usage{
 mediation.test(mv,iv,dv)
}
\arguments{
  \item{mv}{The mediator variable.}
  \item{iv}{The independent variable.}
  \item{dv}{The dependent variable.}
}
\details{
To test whether a mediator carries the influence on an IV to a DV.
}
\value{
Missing values are not allowed.
}
\author{
  B. Wang  \email{bwang@southalabama.edu}
}
\references{
MacKinnon, D. P., & Dwyer, J. H. (1993). Estimating mediated effects in prevention studies. \emph{Evaluation Review}, 17, 144-158.

MacKinnon, D. P., Warsi, G., & Dwyer, J. H. (1995). A simulation study of mediated effect measures. \emph{Multivariate Behavioral Research}, 30, 41-62.

Preacher, K. J., & Hayes, A. F. (2004). SPSS and SAS procedures for estimating indirect effects in simple mediation models. \emph{Behavior Research Methods,Instruments, & Computers}, 36, 717-731.
    
Preacher, K. J., & Hayes, A. F. (2008). asymptotic and resampling strategies for assessing and comparing indirect effects in multiple mediator models. \emph{Behavior Research Methods, Instruments, & Computers}, 40, 879-891.
}
\examples{
mv = rnorm(100)
iv = rnorm(100)
dv = rnorm(100)
mediation.test(mv,iv,dv)
}
\keyword{test}

