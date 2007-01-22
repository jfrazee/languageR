\name{xylowess.fnc}
\alias{xylowess.fnc}
\title{Trellis scatterplot with smoothers}
\description{
  Convenience function for trellis scatterplots with smoothers added.
}
\usage{
xylowess.fnc(fmla = Rating ~ Frequency | Subject, data = weight, 
  span = 2/3, symbolcolor = "darkgrey", 
  linecolor = "blue", xlabel = "", ylabel = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fmla}{A formula.}
  \item{data}{A dataframe.}
  \item{span}{Span for the smoother.}
  \item{symbolcolor}{Color for plot symbols.}
  \item{linecolor}{Color for smoother.}
  \item{xlabel}{Label for horizontal axis.}
  \item{ylabel}{Label for vertical axis.}
}
\value{
  A trellis scatterplot matrix with smoothers is shown on the graphics device.
}
\author{R. H. Baayen}
\seealso{See also \code{\link{xyplot}}.}
\examples{
  data(weightRatings)
  xylowess.fnc(Rating ~ Frequency | Subject, data = weightRatings,
    xlab = "log Frequency", ylab = "Weight Rating")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
