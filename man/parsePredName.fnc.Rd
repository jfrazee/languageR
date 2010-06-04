\name{parsePredName.fnc}
\alias{parsePredName.fnc}
\title{ parse character string specifying restricted cubic spline }
\description{
 parse character string specifying restricted cubic spline into simple
 predictor name and number of knots
}
\usage{
parsePredName.fnc(name)
}
\arguments{
  \item{name}{ character string for predictor, e.g. \code{rcs(X, 3)}}
}
\details{
  not intended for independent use
}
\value{
  a list with components
  \item{baseName}{character string denoting simple predictor name (\code{X})}
  \item{knots}{integer specifying number of knots}
}
%\references{ ~put references to the literature/web site here ~ }
\author{R. H. Baayen}
\note{ 
not intended for independent use
}
\seealso{See Also as \code{\link{plotLMER.fnc}}}
\examples{
  \dontrun{not intended for independent use}
}
\keyword{ regression }
