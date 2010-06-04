\name{degreesOrKnots.fnc}
\alias{degreesOrKnots.fnc}
\title{Extract degree of polynomial or knots for restricted cubic spline}
\description{
  Extract degree of polynomial or knots for restricted cubic spline
  from the predictor name
}
\usage{
degreesOrKnots.fnc(name)
}
\arguments{
  \item{name}{name of predictor, e.g. \code{poly(X, 2, raw = TRUE)}}
}
\details{
  attempts to find degrees or knots if present in input name
}
\value{
  Returns an integer for degrees or knots
}
%\references{None.}
\author{R. H. Baayen}
\note{ 
  not intended for independent use
}
\seealso{See Also as \code{\link{plotLMER.fnc}}}
\examples{
  \dontrun{not intended for independent use}
}
\keyword{ regression }
