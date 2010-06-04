\name{getKnots.fnc}
\alias{getKnots.fnc}
\title{Extracts knots from variable name}
\description{
Extracts knots for predictor specified simply as, e.g., \code{X}
from column names of model@X or model@frame
}
\usage{
getKnots.fnc(colnms, xlb)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{colnms}{ columns of \code{model@X} }
  \item{xlb}{ simple predictor name }
}
\details{
  not intended for independent use
}
\value{
  an integer (number of knots)
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
