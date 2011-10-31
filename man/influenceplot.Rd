% This is a template R document file
% Distributed under GPL 3 or later

\name{influential.plot}
\alias{influential.plot}
\title{Draw plots for the influence measures}
\description{
Draw plots for the influence measures.
}
\usage{
influential.plot(lmobj,type='hadi',ID=FALSE)
}
\arguments{
  \item{lmobj}{An R object by fitting an OLS model to a data set.}
  \item{type}{Plot type.  'hadi': the Hadi's influence Measures; 'potential-residual': potential-residual plot; 'dfits': DFITS plot; 'hat': leverage plot; 'cook': Cook's distance.}
  \item{ID}{Whether to identify points in the plots. Default: FALSE}
}
\value{
  Output the influence measures, including leverage values (Leverage), Hadi's measure (Hadi), Welsch and Kuh Measure (DFIT) and Cook's distance (CookD).  In addition, the standard residuals are also exported.
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
influential.plot(lm0)
influential.plot(lm0,type='hadi')
influential.plot(lm0,type='potential')
influential.plot(lm0,type='leve')
influential.plot(lm0,type='dfit')
influential.plot(lm0,type='cook')
influential.plot(lm0,type='potential',ID=TRUE)


}
\keyword{stats}

