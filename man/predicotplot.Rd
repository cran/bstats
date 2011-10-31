% This is a template R document file
% Distributed under GPL 3 or later

\name{predictor.plot}
\alias{predictor.plot}
\title{Draw plots for predictor impacts on the dependent variable}
\description{
Draw added-variable plot (av) or redidual plus component (rc) plot.
}
\usage{
predictor.plot(lmobj,type='av',ID=FALSE)
}
\arguments{
  \item{lmobj}{An R object by fitting an OLS model to a data set.}
  \item{type}{Plot type.  'av': added variable plot; 'rc': residual plus component plot.}
  \item{ID}{Whether to identify points in the plots. Default: FALSE}
}
\value{
  Missing value not allowed.
}

\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}
\seealso{
 \code{\link{residual.plot}}.
}
\examples{

data(river)
lm0 = lm(Nitrogen~Agr+Forest+Rsdntial+ComIndl, data=river)
predictor.plot(lm0)
predictor.plot(lm0,type='rc')


}
\keyword{stats}

