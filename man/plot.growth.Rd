\name{plot.growth}
\alias{plot.growth}
\title{Plot method for growth objects}
\description{
  This function defines the plot method for growth objects.
}
\usage{
plot.growth(x, w = "all", ...)
}
\arguments{
  \item{x}{A growth object.}
  \item{w}{A character string denoting the name of a specific variable
    to be plotted.}
	\item{...}{other parameters to be passed through to plotting functions.}
}
\value{
  A plot shown on the graphics device.
}
\author{R. H. Baayen}
\seealso{See Also \code{\link{growth.fnc}}.}
\examples{
  data(alice)
  alice.growth = growth.fnc(alice)
  plot(alice.growth)
  plot(alice.growth, w = "Yule")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
